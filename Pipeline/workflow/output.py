from typing import Dict, Optional, Any

from Pipeline.graph.node import Node, NodeTypes
from Pipeline.translations.cwl.cwl import Cwl
from Pipeline.types.data_types import DataType
from Pipeline.tool.tool import ToolOutput, ToolInput


class Output:
    """
        Only catch with output is we infer the type, we don't explicitly define it
    """
    def __init__(self, identifier: str, data_type: Optional[DataType], meta: Any = None,
                 label: str=None, doc: str=None):
        self._identifier: str = identifier
        self.label = label
        self.doc = doc

        self.data_type: DataType = data_type
        self.meta = meta

    def id(self) -> str:
        return self._identifier

    def cwl(self, output_source: str):
        d = {
            **self.data_type.cwl(),
            Cwl.WORKFLOW.OUTPUT.kID: self._identifier,
            Cwl.WORKFLOW.OUTPUT.kOUTPUT_SOURCE: output_source
        }

        if self.label:
            d[Cwl.WORKFLOW.kLABEL] = self.label
        if self.doc:
            d[Cwl.WORKFLOW.kDOC] = self.doc

        return d


class OutputNode(Node):

    def __init__(self, output: Output):
        self.output = output
        super().__init__(NodeTypes.OUTPUT, output.id())

    def inputs(self) -> Dict[str, ToolInput]:
        return {self.output._identifier: ToolInput(self.output.meta, self.output.data_type)}

    def outputs(self) -> Dict[str, ToolOutput]:
        return {}

    def cwl(self):
        source = next(iter(self.connection_map.values()))
        return self.output.cwl(source[0])
