from typing import Dict

from pipeline_definition.types.data_types import DataType
from pipeline_definition.types.step import ToolOutput
from pipeline_definition.graph.node import Node, NodeType


class Input:
    def __init__(self, label: str, data_type: DataType):
        self.label: str = label
        self.data_type: DataType = data_type

    def id(self):
        return self.label


class InputNode(Node):

    def __init__(self, inp: Input):
        super().__init__(NodeType.INPUT, inp.id())
        self.input: Input = inp

    def outputs(self) -> Dict[str, ToolOutput]:
        # Program will just grab first value anyway
        return {"output": ToolOutput(self.input.label, self.input.data_type)}

    def inputs(self):
        return None

    def translate(self, mapped_inputs):
        ind = dict()
        for inp in self.__workflowInputSet:
            ind.update(inp.translate_for_workflow())
        return {'inputs': ind}

    @staticmethod
    def input_step_tag_name():
        return 'input'

