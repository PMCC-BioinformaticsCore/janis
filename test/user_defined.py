from pipeline_definition.types.step import StepFactory
from pipeline_definition.types.step import Step


class UserDefinedStepFactory(StepFactory):
    @classmethod
    def type(cls):
        return 'userdefined-step'

    @classmethod
    def label(cls):
        return 'userdefined-step'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def schema(cls):
        return {
            'schema': {
                'aligner': {
                    'type': 'string'
                }
            },
            'nullable': True
        }

    @classmethod
    def build(cls, meta):
        step = UseerDefinedStep( meta )
        return step

class UseerDefinedStep(Step):
    pass







