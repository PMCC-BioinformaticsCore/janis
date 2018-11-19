from typing import Any, Dict

from pipeline_definition.types.input_type import InputFactory, InputType
from pipeline_definition.types.input_type import Input
from pipeline_definition.utils.logger import Logger

class_file = InputType('class', label='a java .class file', description='A compiled java binary class file')


class ClassFile(Input):
    def translate_for_workflow(self) -> dict:
        raise Exception('Not yet implemented')

    def translate_for_input(self):
        fd = {'class': 'File', 'path': self._path}
        return {self.id(): fd}

    def resolve(self):
        pass

    def __init__(self, label: str, meta: Dict[str, Any]):
        super().__init__(label, meta)
        self._resolved = True
        self._path = [self.meta()['path']]
        self._resolved = False

    def identify(self):
        Logger.log(f" {self.id()} :: ClassFile: {self._path}")

    def datum_type(self):
        return self.type()

    def is_subtype_of(self, other):
        return False


class ClassFileFactory(InputFactory):
    @classmethod
    def type(cls) -> InputType:
        return class_file

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
    def build(cls, label: str, meta: Dict[str, Any]) -> ClassFile:
        return ClassFile(label, meta)
