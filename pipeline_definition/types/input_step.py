from typing import Dict

from pipeline_definition.types.input_type import InputType, get_input_type
from pipeline_definition.types.step_type import Step


class InputStep(Step):
  def translate(self, mapped_inputs):
    ind = dict()
    for inp in self.__workflowInputSet:
      ind.update(inp.translate_for_workflow())
    return {'inputs': ind}

  def __init__(self, workflow_input_set):
    super().__init__({
      'inputs': {
        'WorkflowInputStep': {
        },
        'tag': InputStep.input_step_tag_name()
      }})
    self.__workflowInputSet = workflow_input_set

  def provides(self) -> Dict[str, InputType]:
    if not self.__workflowInputSet:
      return {}

    outputs = dict()
    for inp in self.__workflowInputSet:
      outputs[inp.id()] = inp.type().type_name()

    return outputs

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
