from pipeline_definition.types.input_type import InputFactory
from pipeline_definition.types.input_type import Input

class PairedReadFactory(InputFactory):
    @classmethod
    def type(cls):
        return 'SequenceReadArchivePaired'

    @classmethod
    def label(cls):
        return 'paired read files'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def describe(cls):
        return {
            'schema' : {
                'forward-pattern': {'type': 'string'},
                'backward-pattern': {'type': 'string'}
            }
        }

    @classmethod
    def build(cls, meta):
        input = PairedReadInput( meta )
        return input


class PairedReadInput(Input):
    pass

