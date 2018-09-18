from pipeline_definition.types.step_type import Step


class InputStep(Step):
  def translate(self):
    pass

  def __init__(self, workflow_input_set):
    super().__init__({
      'input-step': {
        'WorkflowInputStep': {
        },
        'tag': InputStep.input_step_tag_name()
      }})
    self.__workflowInputSet = workflow_input_set

  def provides(self):
    # return self.__workflowInputSet
    if not self.__workflowInputSet:
      return None

    outputs = list()
    for input in self.__workflowInputSet:
      o = {
        Step.STR_ID: input.id(),
        Step.STR_TYPE: input.type()
      }
      if not outputs:
        outputs = list()

      outputs.append(o)

    return outputs

  def requires(self):
    return None

  def inputs(self):
    return self.__workflowInputSet

  def type(self):
    return 'input'

  @staticmethod
  def input_step_tag_name():
    return 'input'
