from typing import Dict, List

from pipeline_definition.types.input_type import InputType, Input
from pipeline_definition.types.step_type import Step
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

    # super().__init__({
    #   'inputs': {
    #     'WorkflowInputStep': {
    #     },
    #     'tag': InputNode.input_step_tag_name()
    #   }})

  def provides(self) -> InputType:
    return self.input.type()
    # return { inp.id(): inp.type() for inp in self.__workflowInputSet}
    # for inp in self.__workflowInputSet:
    #   outputs[inp.id()] = get_input_type(inp.type().type_name())
    # return outputs

  def requires(self):
    return None

  def inputs(self):
    return self.__workflowInputSet

  def type(self):
    return 'input'

  def cores(self):
    return None

  def ram(self):
    return None

  @staticmethod
  def input_step_tag_name():
    return 'input'
