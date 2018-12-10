from typing import Dict, Optional, Any

from Pipeline.graph.node import Node, NodeTypes
from Pipeline.types.data_types import DataType
from Pipeline.tool.tool import ToolOutput, ToolInput


class Output:
    """
        Only catch with output is we infer the type, we don't explicitly define it
    """
    def __init__(self, label: str, data_type: Optional[DataType], meta: Any = None):
        self.label: str = label
        self.data_type: DataType = data_type
        self.meta = meta

    def id(self) -> str:
        return self.label


class OutputNode(Node):

    def __init__(self, output: Output):
        self.output = output
        super().__init__(NodeTypes.OUTPUT, output.label)

    def inputs(self) -> Dict[str, ToolInput]:
        return {self.output.label: ToolInput(self.output.meta, self.output.data_type)}

    def outputs(self) -> Dict[str, ToolOutput]:
        return {}

    def cwl(self):
        return {
            **self.output.data_type.cwl(),
            "outputSource": next(iter(self.connection_map.values()))[0]
        }
