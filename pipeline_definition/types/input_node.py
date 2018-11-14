from typing import Dict, List

from pipeline_definition.types.input_type import InputType, Input
from pipeline_definition.types.step_type import Step, StepOutput
from pipeline_definition.graph.node import Node, NodeType


class InputNode(Node):
    def translate(self, mapped_inputs):
        ind = dict()
        for inp in self.__workflowInputSet:
            ind.update(inp.translate_for_workflow())
        return {'inputs': ind}

    def __init__(self, inp: Input):
        super().__init__(NodeType.INPUT, inp.id())
        self.input: Input = inp

    def outputs(self) -> Dict[str, StepOutput]:
        # Program will just grab first value anyway
        return {"output": StepOutput(self.input.id(), self.input.type())}

    def inputs(self):
        return None

    @staticmethod
    def input_step_tag_name():
        return 'input'
