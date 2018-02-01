from pipeline_definition.types.input_type import InputFactory
from pipeline_definition.types.input_type import Input

class BAMFactory(InputFactory):
    @classmethod
    def type(cls):
        return 'BAM'

    @classmethod
    def label(cls):
        return 'BAM file'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def describe(cls):
        return {
            'schema' : {
                'path': {'type': 'string'},
                'label': {'type': 'string'}
            },
            'nullable': True

        }

    @classmethod
    def build(cls, dict):
        input = BAMInput( dict )
        return input


class BAMInput(Input):
    def __init__(self, dict):
        super().__init__( dict )
        self.path = None

        if self.meta is not None:
            self.path = self.meta.get("path")

    def identify(self):
        super().identify()
        print("Path:", self.path)




