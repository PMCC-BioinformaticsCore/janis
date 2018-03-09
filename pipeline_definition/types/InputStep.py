
from pipeline_definition.types.step_type import Step


class InputStep(Step):
    def __init__(self, workflowInputSet):
        super().__init__({
            'input-step' : {
                'WorkflowInputStep' : {
                },
                'tag' : InputStep.inputSteptagName()
            }})
        self.__workflowInputSet = workflowInputSet

    def provides(self):
        if not self.__workflowInputSet:
            return None

        outputs = None
        for input in self.__workflowInputSet:
            o = {
                Step.STR_ID : input.id(),
                Step.STR_TYPE : input.type()
            }
            if not outputs:
                outputs = list()

            outputs.append(o)


        return outputs

    def requires(self):
        return None

    @staticmethod
    def inputSteptagName():
        return 'input'