from typing import Dict, Optional, Any

from janis.utils.validators import Validators

from janis.graph.node import Node, NodeTypes
from janis.types.data_types import DataType
from janis.tool.tool import ToolOutput, ToolInput


class Output:
    """
        Only catch with output is we infer the type, we don't explicitly define it
    """

    illegal_keywords = ["output"]

    def __init__(
        self, identifier: str, meta: Any = None, label: str = None, doc: str = None
    ):
        self._identifier: str = identifier

        if not Validators.validate_identifier(identifier):
            raise Exception(
                f"The output identifier '{identifier}' was not validated by '{Validators.identifier_regex}'"
                f" (must start with letters, and then only contain letters, numbers and an underscore)"
            )

        if identifier in Output.illegal_keywords:
            raise Exception(
                f"The output identifier '{identifier}' is a reserved keyword "
                f"({', '.join(Output.illegal_keywords)})"
            )

        self.label = label
        self.doc = doc

        self.data_type: Optional[DataType] = None
        self.meta = meta

    def id(self) -> str:
        return self._identifier


class OutputNode(Node):
    def __init__(self, output: Output):
        self.output = output
        super().__init__(NodeTypes.OUTPUT, output.id())

    def inputs(self) -> Dict[str, ToolInput]:
        if not self.output.data_type:
            raise Exception(
                f"Output '{self.output.id()}' did not resolve with a DataType, "
                f"this likely means it wasn't connected to step or an internal error has occurred."
            )

        return {self.output.id(): ToolInput("outp", self.output.data_type)}

    def outputs(self) -> Dict[str, ToolOutput]:
        return {}
