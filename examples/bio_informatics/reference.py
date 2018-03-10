from pipeline_definition.types.input_type import InputFactory
from pipeline_definition.types.input_type import Input

class ReferenceFactory(InputFactory):
    @classmethod
    def type(cls):
        return 'REFERENCE'

    @classmethod
    def label(cls):
        return 'A reference genome'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def describe(cls):
        return {
            'schema': {
                'path': {'type': 'string'},
                'label': {'type': 'string'}
            },
            'nullable': True

        }

    @classmethod
    def build(cls, dict):
        input = ReferenceInput( dict )
        return input


class ReferenceInput(Input):
    def __init__(self, dict):
        super().__init__( dict )
        self.path = None

        if self.meta is not None:
            self.path = self.meta().get("path")

    def identify(self):
        super().identify()
        print("Path:", self.path)

    def datum_type(self):
        return self.type()

    def is_subtype_of(self, other):
        return False



