
from pipeline_definition.types.step_type import Step


class InputStep(Step):
    def __init__(self, workflowInputSet):
        super().__init__({
            'root' : {
                'WorkflowInputStep' : {
                },
                'tag' : 'root'
            }})
        self.__workflowInputSet = workflowInputSet

    def provides(self):
        pass

    def requires(self):
        return None
