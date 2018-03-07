from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step

class CallFactory(StepFactory):
    @classmethod
    def type(cls):
        return 'call'

    @classmethod
    def label(cls):
        return 'call'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def describe(cls):
        return {
            'schema': {
                'caller': {
                    'type': 'string'
                }
            },
            'nullable': True
        }

    @classmethod
    def build(cls, meta):
        step = CallStep( meta )
        return step

class CallStep(Step):

    def provides(self):
        pass

    def requires(self):
        pass
