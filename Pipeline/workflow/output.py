from typing import Dict

from Pipeline.graph.node import Node, NodeTypes
from Pipeline.types.data_types import DataType
from Pipeline.tool.tool import ToolOutput, ToolInput


class Output:
    """
        Only catch with output is we infer the type, we don't explicitly define it
    """
    def __init__(self, label: str, data_type: DataType):
        self.label: str = label
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
        super().__init__(NodeTypes.OUTPUT, output.label)

    def inputs(self) -> Dict[str, ToolInput]:
        return {self.output.label: ToolInput(self.output.source, self.output.data_type)}

    def outputs(self) -> Dict[str, ToolOutput]:
        return {}
