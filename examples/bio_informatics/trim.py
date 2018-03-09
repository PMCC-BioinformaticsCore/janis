from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step

class TrimFactory(StepFactory):
    @classmethod
    def type(cls):
        return 'trim'

    @classmethod
    def label(cls):
        return 'trim'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def describe(cls):
        return {
            'schema': {
                'trimmer': {
                    'type': 'string',
                    'allowed': ['cutadapt', 'trimmomatic'],
                    'default': 'trimmomatic'
                }
            },
            'nullable': True
        }

    @classmethod
    def build(cls, meta):
        step = TrimStep( meta )
        return step

class TrimStep(Step):

    def provides(self):
        return {
            "type": "SequenceReadArchive"
        }


    def requires(self):
        return {
            "type": "SequenceReadArchive"
        }

