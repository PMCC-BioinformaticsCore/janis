from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step

class JointCallFactory(StepFactory):
    @classmethod
    def type(cls):
        return 'joint_call'

    @classmethod
    def label(cls):
        return 'joint_call'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def describe(cls):
        return {
            'schema': {
                'normal_tag': {
                    'type': 'string'
                },
                'tumour_tag': {
                    'type': 'string'
                }
            },
            'nullable': True
        }

    @classmethod
    def build(cls, meta):
        step = JointCallStep( meta )
        return step

class JointCallStep(Step):

    def provides(self):
        pass

    def requires(self):
        pass
