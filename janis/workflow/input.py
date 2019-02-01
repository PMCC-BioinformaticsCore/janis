from typing import Dict, Optional, List

from janis.utils.logger import Logger, LogLevel
# import cwlgen.cwlgen as cwl
from janis.types.data_types import DataType, NativeTypes
from janis.workflow.step import ToolOutput
from janis.graph.node import Node, NodeTypes


class Input:
    def __init__(self, identifier: str, data_type: DataType, value_from_input=None,
                 label: str=None, doc: str=None):
        Logger.log(f"Creating input '{identifier}' with value: '{value_from_input}'")
        self._identifier: str = identifier

        self.label = label
        self.doc = doc

        self.value = value_from_input

        self.data_type: DataType = data_type

    def id(self):
        return self._identifier

    def cwl_input(self):
        return self.data_type.cwl_input(self.value)

    def wdl_input(self):
        return self.value


class InputNode(Node):

    def __init__(self, inp: Input):
        Node.__init__(self, NodeTypes.INPUT, inp.id())
        self.input: Input = inp

    def outputs(self) -> Dict[str, ToolOutput]:
        # Program will just grab first value anyway
        return {"": ToolOutput(self.input._identifier, self.input.data_type)}

    def inputs(self):
        return None
