from pipeline_definition.types.step_type import Step


class InputStep(Step):
  def __init__(self, workflow_input_set):
    super().__init__({
      'input-step': {
        'WorkflowInputStep': {
        },
        'tag': InputStep.input_step_tag_name()
      }})
    self.__workflowInputSet = workflow_input_set

  def provides(self):
    if not self.__workflowInputSet:
      return None

    outputs = None
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

  @staticmethod
  def input_step_tag_name():
    return 'input'
