from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step

class DedupFactory(StepFactory):
    @classmethod
    def type(cls):
        return 'dedup'

    @classmethod
    def label(cls):
        return 'dedup'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def describe(cls):
        return {
            'schema': {
                'dedup': {
                    'type': 'string'
                }
            },
            'nullable': True
        }

    @classmethod
    def build(cls, meta):
        step = DedupStep( meta )
        return step

class DedupStep(Step):

    def provides(self):
        pass

    def requires(self):
        pass
