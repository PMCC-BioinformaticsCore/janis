
from pipeline_definition.types.step_type import Step


class InputStep(Step):
    def __init__(self, workflowInputSet):
        super().__init__({
            'input-step' : {
                'WorkflowInputStep' : {
                },
                'tag' : InputStep.tagConvention()
            }})
        self.__workflowInputSet = workflowInputSet

    def provides(self):

        if not self.__workflowInputSet:
            return None








        return self.__workflowInputSet

    def requires(self):
        return None

    @staticmethod
    def tagConvention():
        return 'default-branch'