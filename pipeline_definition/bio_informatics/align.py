from pipeline_definition.types.step_type import StepFactory

class AlignFactory(StepFactory):
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
        print(">>>>>>>>>>>>>>>> ", meta )
        return None

    @classmethod
    def type(cls):
        return 'align'

    @classmethod
    def label(cls):
        return 'Aligner'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def emit(cls):
        return "Fcatory says: " + cls.__name__

