from typing import Dict

from pipeline_definition.graph.node import Node, NodeType
from pipeline_definition.types.data_types import DataType
from pipeline_definition.types.tool import ToolOutput, ToolInput


class Output:
    """
        Only catch with output is we infer the type, we don't explicitly define it
    """
    def __init__(self, label: str, source: str, data_type: DataType):
        self.label: str = label
        self.source: str = source
        self.data_type: DataType = data_type

    def id(self):
        return self.label

    def cwl(self):
        return {
            **self.data_type.cwl(),
            "outputSource": self.source
        }


class OutputNode(Node):

    def __init__(self, output: Output):
        self.output = output
        super().__init__(NodeType.OUTPUT, output.label)

    def inputs(self) -> Dict[str, ToolInput]:
        return {self.output.label: ToolInput(self.output.source, self.output.data_type)}

    def outputs(self) -> Dict[str, ToolOutput]:
        return {}
