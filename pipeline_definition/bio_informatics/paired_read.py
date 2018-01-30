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
                'label' : {'type': 'string'},
                'forward-pattern': {'type': 'string', 'required': True },
                'backward-pattern': {'type': 'string'}
            }
        }

    @classmethod
    def build(cls, dict):
        input = PairedReadInput( dict )
        return input


class PairedReadInput(Input):
    def __init__(self, dict):
        super().__init__( dict )
        self.forwardPattern = self.meta["forward-pattern"]
        self.backwardPattern = self.meta["backward-pattern"]

    def identify(self):
        super().identify()
        print("Forward Pattern:", self.forwardPattern)
        print("Backward Pattern:", self.backwardPattern)
