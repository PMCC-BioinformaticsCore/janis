from typing import Optional, Dict, Any


from janis.graph.node import Node, NodeTypes
from janis.tool.tool import ToolOutput, ToolInput
from janis.types.common_data_types import Array
from janis.utils import first_value
from janis.utils.logger import Logger


def full_lbl(node: Node, tag: Optional[str]) -> str:
    if tag is None:
        return node.id()
    return f"{node.id()}/{tag}"


def full_dot(node: Node, tag: Optional[str]) -> str:
    if tag is None:
        return node.id()
    return f"{node.id()}.{tag}"


class Edge:
    def __init__(
        self, start: Node, stag: Optional[str], finish: Node, ftag: Optional[str]
    ):
        Logger.log(
            f"Creating edge: ({NodeTypes.to_str(start.node_type)}) '{start.id()}.{stag}' → "
            f"({NodeTypes.to_str(finish.node_type)}) '{finish.id()}.{ftag}'"
        )

        self.start: Node = start
        self.stag: Optional[str] = stag
        self.finish: Node = finish
        self.ftag: Optional[str] = ftag
        self.compatible_types: Optional[bool] = None
        self.scatter = None

        self.validate_tags()
        self.check_types()

    def source(self):
        return full_lbl(self.start, self.stag)

    def dotted_source(self):
        return full_dot(self.start, self.stag)

    def validate_tags(self):
        if (
            self.start.node_type == NodeTypes.TASK
            and self.stag not in self.start.outputs()
        ):
            raise Exception(
                f"Could not find the tag '{self.stag}' in the inputs of '{self.start.id()}'"
            )
        if (
            self.finish.node_type == NodeTypes.TASK
            and self.ftag not in self.finish.inputs()
        ):
            raise Exception(
                f"Could not find the tag '{self.ftag}' in the outputs of '{self.finish.id()}': {list(self.finish.inputs().keys())}"
            )

    def check_types(self):
        from janis.workflow.input import InputNode

        stype: ToolOutput = self.start.outputs()[
            self.stag
        ] if self.stag is not None else first_value(self.start.outputs())
        ftype: ToolInput = self.finish.inputs()[
            self.ftag
        ] if self.ftag is not None else first_value(self.finish.inputs())

        self.compatible_types = False
        self.scatter = False

        start_is_scattered = any(
            e.has_scatter() for e in self.start.connection_map.values()
        )

        source_has_default = (
            isinstance(self.start, InputNode) and self.start.input.default is not None
        )
        if not start_is_scattered and ftype.input_type.can_receive_from(
            stype.output_type, source_has_default=source_has_default
        ):
            self.compatible_types = True
            self.scatter = False
            return

        if isinstance(ftype.input_type, Array):
            # This might be a merge

            # check if s has a scatter step, then we sweet
            if start_is_scattered and self.finish.node_type != NodeTypes.OUTPUT:
                Logger.log(
                    f"This edge merges the inputs from '{full_dot(self.start, self.stag)}' for "
                    f"'{full_dot(self.finish, self.ftag)}'"
                )
                self.compatible_types = ftype.input_type.subtype().can_receive_from(
                    stype.output_type
                )
            else:
                self.compatible_types = ftype.input_type.can_receive_from(
                    stype.output_type
                )

            # Check if its array -> single (scatter)
            if not self.compatible_types and ftype.input_type.subtype().can_receive_from(
                stype.output_type
            ):
                self.compatible_types = True
                s = full_dot(self.start, self.stag)
                f = full_dot(self.finish, self.ftag)
                Logger.log(f"Safely provided '{s}' to the array input '{f}'")

        if isinstance(stype.output_type, Array):
            # Potential scattering event - non-scatter if s is compatible, scatter if s.subtype is compatible

            # Scattering event if ftype.input_type.canReceiveFrom(stype.output_type.subtype)
            if ftype.input_type.can_receive_from(stype.output_type.subtype()):
                Logger.info(
                    f"Scatter the connection between '{full_dot(self.start, self.stag)}' → "
                    f"'{full_dot(self.finish, self.ftag)}' ({stype.output_type.id()} → {ftype.input_type.id()})"
                )
                self.compatible_types = True
                self.scatter = True
                return

        if start_is_scattered:
            # Potential rescattering event

            # Scattering event if ftype.input_type.canReceiveFrom(stype.output_type.subtype)
            if ftype.input_type.can_receive_from(stype.output_type):
                Logger.info(
                    f"Rescatter the connection between '{full_dot(self.start, self.stag)}' → "
                    f"'{full_dot(self.finish, self.ftag)}' ({stype.output_type.id()} → {ftype.input_type.id()})"
                )
                self.compatible_types = True
                self.scatter = True
                return

        if not self.compatible_types:
            s = full_dot(self.start, self.stag)
            f = full_dot(self.finish, self.ftag)
            Logger.critical(
                f"Mismatch of types when joining '{s}' to '{f}' "
                f"({stype.output_type.id()} -/→ {ftype.input_type.id()})"
            )
            Logger.log(f"No action taken to correct type-mismatch of '{s}' to {f}'")


class StepInput:
    def __init__(self, finish: Node, finish_tag: str):

        self.finish: Node = finish
        self.ftag: Optional[str] = finish_tag

        self.default = None
        self.multiple_inputs = False

        self.source_map: Dict[str, Edge] = {}

    def add_source(self, start: Node, stag: Optional[str]):
        """
        Add a connection
        :param start:
        :param stag:
        :return:
        """
        finish_type = (
            self.finish.inputs()[self.ftag]
            if self.ftag is not None
            else first_value(self.finish.inputs())
        ).input_type

        if len(self.source_map) == 1 and start.id() not in self.source_map:
            self.multiple_inputs = True

            if not isinstance(finish_type, Array):
                Logger.warn(
                    f"Adding multiple inputs to '{self.finish.id()}' and '{finish_type.id()}' is not an array"
                )

        e = Edge(start, stag, self.finish, self.ftag)
        self.source_map[start.id()] = e
        return e

    def has_scatter(self):
        return any(e.scatter for e in self.source_map.values())

    def set_default(self, default: Any):
        Logger.log(
            f"Setting the default of '{self.finish.id()}.{self.ftag}' to be '{str(default)}'"
        )
        self.default = default

    def source(self):
        n = len(self.source_map)
        if n == 0:
            return None
        elif n == 1:
            return first_value(self.source_map)
        else:
            return list(self.source_map.values())

    def dotted_source(self):
        n = len(self.source_map)
        if n == 0:
            return None
        elif n == 1:
            return first_value(self.source_map).dotted_source()
        else:
            return [e.dotted_source() for e in self.source_map.values()]

    def slashed_source(self):
        n = len(self.source_map)
        if n == 0:
            return None
        elif n == 1:
            return first_value(self.source_map).source()
        else:
            return [e.source() for e in self.source_map.values()]
