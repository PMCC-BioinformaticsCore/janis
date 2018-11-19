from pipeline_definition.types.input_type import InputFactory, InputType
from pipeline_definition.types.input_type import Input
from typing import Dict

generic_string = InputType('string', label='a generic string')


class String(Input):
    def translate_for_workflow(self) -> dict:
        raise Exception('Not yet implemented')

    def translate_for_input(self):
        return {self.id(): self._value}

    def resolve(self):
        pass

    def __init__(self, label: str, meta: Dict):
        # meta will actually be a string
        super().__init__(label, meta)
        self._value = str(meta)

    def identify(self):
        super().identify()

    def datum_type(self):
        return self.type()

    def is_subtype_of(self, other):
        return False


class StringFactory(InputFactory):
    @classmethod
    def type(cls) -> InputType:
        return generic_string

    @classmethod
    def schema(cls):
        return {
            'schema': {
                'path': {'type': 'string'},
                'label': {'type': 'string'}
            },
            'nullable': True
        }

    @classmethod
    def build(cls, label: str, meta: Dict) -> String:
        return String(label, meta)
