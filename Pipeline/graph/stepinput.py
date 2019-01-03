from typing import Tuple, Optional, Dict, List, Any

from Pipeline.utils.logger import Logger
from Pipeline.tool.commandtool import ToolOutput, ToolInput
from Pipeline.types.common_data_types import Array
from Pipeline.graph.node import Node, NodeTypes
from Pipeline.workflow.step import StepNode


def first_value(d: Dict):
    return next(iter(d.values()))


def full_lbl(node: Node, tag: Optional[str]) -> str:
    if tag is None:
        return node.id()
    return "/".join([node.id(), tag])


class Edge:
    def __init__(self, start: Node, stag: Optional[str], finish: Node, ftag: Optional[str]):
        Logger.log(f"Creating edge: ({NodeTypes.to_str(start.node_type)}) '{start.id()}.{stag}' → "
                   f"({NodeTypes.to_str(finish.node_type)}) '{finish.id()}.{ftag}'")

        self.start: Node = start
        self.stag: Optional[str] = stag
        self.finish: Node = finish
        self.ftag: Optional[str] = ftag
        self.compatible_types: bool = None
        self.scatter = None

        self.validate_tags()
        self.check_types()

    def source(self):
        return full_lbl(self.start, self.stag)

    def validate_tags(self):
        if self.start.node_type == NodeTypes.TASK and self.stag not in self.start.outputs():
            raise Exception(f"Could not find the tag '{self.stag}' in the inputs of '{self.start.id()}'")
        if self.finish.node_type == NodeTypes.TASK and self.ftag not in self.finish.inputs():
            raise Exception(f"Could not find the tag '{self.ftag}' in the outputs of '{self.finish.id()}': {list(self.finish.inputs().keys())}")

    def check_types(self):
        stype: ToolOutput = self.start.outputs()[self.stag] \
            if self.stag is not None else first_value(self.start.outputs())
        ftype: ToolInput = self.finish.inputs()[self.ftag] \
            if self.ftag is not None else first_value(self.finish.inputs())

        s_is_array = isinstance(stype.output_type, Array)
        f_is_array = isinstance(ftype.input_type, Array)

        if s_is_array and not f_is_array:
            Logger.log(f"Potential scattering event between '{self.start.id()}' ({stype.output_type.id()}) "
                       f"and '{self.finish.id()}' ({ftype.input_type.id()})")
            # this second check means that the mypy (and linter) is happy with .subtype()
            self.compatible_types = isinstance(stype.output_type, Array) and \
                                    ftype.input_type.can_receive_from(stype.output_type.subtype())
        elif f_is_array and not s_is_array:
            # check if s has a scatter step, then we sweet
            start_is_scattered = any(e.scatter for e in self.start.connection_map.values())
            if start_is_scattered and isinstance(ftype.input_type, Array):
                self.compatible_types = ftype.input_type.subtype().can_receive_from(stype.output_type)
            else:
                self.compatible_types = ftype.input_type.can_receive_from(stype.output_type)
        else:
            self.compatible_types = ftype.input_type.can_receive_from(stype.output_type)

        if not self.compatible_types:
            s = full_lbl(self.start, self.stag)
            f = full_lbl(self.finish, self.ftag)
            Logger.critical(f"Mismatch of types when joining '{s}' to '{f}' "
                            f"({stype.output_type.id()} -/→ {ftype.input_type.id()})")
            Logger.log(f"No action taken to correct type-mismatch of '{s}' to {f}'")

        self.scatter: bool = s_is_array and not f_is_array and self.compatible_types


class StepInput:
    # todo: Make Edge support multiple source_map

    def __init__(self, finish: Node, finish_tag: str):

        self.finish: Node = finish
        self.ftag: Optional[str] = finish_tag

        self.default = None
        self.scatter = False
        self.multiple_inputs = False

        self.source_map: Dict[str, Edge] = {}

    def add_source(self, start: Node, stag: Optional[str]):
        """
        Add a connection
        :param start:
        :param stag:
        :return:
        """
        finish_type = (self.finish.inputs()[self.ftag]
                       if self.ftag is not None else
                       first_value(self.finish.inputs())).input_type

        if len(self.source_map) == 1 and start.id() not in self.source_map:
            self.multiple_inputs = True

            if not isinstance(finish_type, Array):
                Logger.warn(f"Adding multiple inputs to '{self.finish.id()}' and '{finish_type.id()}' is not an array")

        e = Edge(start, stag, self.finish, self.ftag)
        self.source_map[start.id()] = e
        return e

    def set_default(self, default: Any):
        Logger.log(f"Setting the default of '{self.finish.id()}.{self.ftag}' to be '{str(default)}'")
        self.default = default

    def source(self):
        n = len(self.source_map)
        if n == 0:
            return None
        elif n == 1:
            return first_value(self.source_map).source()
        else:
            return [e.source() for e in self.source_map.values()]
