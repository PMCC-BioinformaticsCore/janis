from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step

class AlignFactory(StepFactory):
    @classmethod
    def type(cls):
        return 'align'

    @classmethod
    def label(cls):
        return 'align'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def describe(cls):
        return {
            'schema': {
                'aligner': {
                    'type': 'string',
                    'allowed': ['bowtie', 'bwa'],
                    'default': 'bowtie'
                }
            },
            'nullable': True
        }

    @classmethod
    def build(cls, meta):
        step = AlignStep( meta )
        return step

class AlignStep(Step):

    def provides(self):
        return [
            {
                Step.STR_ID: "alignedbamfile",
                Step.STR_TYPE: "BAM"
            }
        ]

    def requires(self):
        return [
            {
                Step.STR_ID: "read",
                Step.STR_TYPE: "SequenceReadArchivePaired"
            }
        ]
