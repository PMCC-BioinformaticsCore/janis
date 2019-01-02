from typing import Tuple, Optional, Dict, List

from Pipeline.utils.logger import Logger
from Pipeline.tool.commandtool import ToolOutput, ToolInput
from Pipeline.types.common_data_types import Array
from Pipeline.graph.node import Node, NodeTypes
from Pipeline.workflow.step import StepNode

TaggedNode = Tuple[Node, Optional[List[str]]]


def first_value(d: Dict):
    return next(iter(d.values()))


def full_lbl(node: Node, tag: Optional[str]) -> str:
    if tag is None:
        return node.id()
    return "/".join([node.id(), tag])


class Edge:

    def __init__(self, start: Optional[TaggedNode], finish: TaggedNode, default=None):
        if start:
            self.start: Node = start[0]
            self.stag: Optional[str] = start[1][1] if len(start[1]) > 1 else None
        else:
            self.start, self.stag = None, None

        # The way we'll store an edge, we won't require this
        self.finish: Node = finish[0]
        self.ftag: Optional[str] = finish[1][1] if len(finish[1]) > 1 else None

        self.default = default
        self.scatter = False
        self.correct_type = False

        if self.start:
            Logger.log(f"Creating edge: ({NodeTypes.to_str(self.start.node_type)}) '{self.start.id()}.{self.stag}' → "
                       f"({NodeTypes.to_str(self.finish.node_type)}) '{self.finish.id()}.{self.ftag}'")
        else:
            Logger.log(f"Created empty entry point to '{self.finish.id()}.{self.ftag}' "
                       f"({NodeTypes.to_str(self.finish.node_type)})")

        if self.default:
            Logger.log(f"Added default value: '{default}' to previously created edge")

        if self.stag is not None and self.stag not in self.start.outputs():
            keys = ", ".join(self.start.outputs().keys())
            tool = ""
            if isinstance(self.start, StepNode):
                tool = f"(Tool: {self.finish.step.tool().id()}) "
            raise Exception(f"Could not find the tag '{self.stag}' as an input of '{self.start.id()}', "
                            f"{tool}expected one of: {keys}")
        if self.ftag is not None and self.ftag not in self.finish.inputs():
            keys = ", ".join(self.finish.inputs().keys())
            tool = ""
            if isinstance(self.finish, StepNode):
                tool = f"(Tool: {self.finish.step.tool().id()}) "
            raise Exception(f"Could not find the tag '{self.stag}' as an output of '{self.finish.id()}', "
                            f"{tool}expected one of: {keys}")

        if self.start and self.finish:
            stype: ToolOutput = self.start.outputs()[self.stag] \
                if self.stag is not None else first_value(self.start.outputs())
            ftype: ToolInput = self.finish.inputs()[self.ftag] \
                if self.ftag is not None else first_value(self.finish.inputs())

            s_is_array = isinstance(stype.output_type, Array)
            f_is_array = isinstance(ftype.input_type, Array)

            self.correct_type = None
            if s_is_array and not f_is_array:
                Logger.log(f"Potential scattering event between '{self.start.id()}' ({stype.output_type.id()}) "
                           f"and '{self.finish.id()}' ({ftype.input_type.id()})")
                # this second check means that the mypy (and linter) is happy with .subtype()
                self.correct_type = isinstance(stype.output_type, Array) and \
                              ftype.input_type.can_receive_from(stype.output_type.subtype())
            elif f_is_array and not s_is_array:
                # check if s has a scatter step, then we sweet
                start_is_scattered = any(e.scatter for e in self.start.connection_map.values())
                if start_is_scattered and isinstance(ftype.input_type, Array):
                    self.correct_type = ftype.input_type.subtype().can_receive_from(stype.output_type)
                else:
                    self.correct_type = ftype.input_type.can_receive_from(stype.output_type)
            else:
                self.correct_type = ftype.input_type.can_receive_from(stype.output_type)

            if not self.correct_type:
                s = full_lbl(self.start, self.stag)
                f = full_lbl(self.finish, self.ftag)
                Logger.critical(f"Mismatch of types when joining '{s}' to '{f}' "
                                f"({stype.output_type.id()} -/→ {ftype.input_type.id()})")
                Logger.log(f"No action taken to correct type-mismatch of '{s}' to {f}'")

            self.scatter: bool = s_is_array and not f_is_array and self.correct_type
        else:
            Logger.log("Didn't have start or finish: skipping type check for this edge")

    def source(self):
        if self.start:
            return full_lbl(self.start, self.stag)
        return None
