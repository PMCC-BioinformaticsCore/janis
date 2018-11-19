# from pipeline_definition.types.input_type import InputFactory, InputType
from pipeline_definition.types.input_type import Input
from typing import Dict

# generic_file = InputType('file', label='an untyped file')


# class GenericFile(Input):
#   def translate_for_workflow(self) -> dict:
#     raise Exception('Not yet implemented')
#
#   def translate_for_input(self):
#     fd = {'class': 'File', 'path': self._path}
#     return {self.id(): fd}
#
#   def resolve(self):
#     pass
#
#   def __init__(self, label: str, meta: Dict, debug=False):
#     super().__init__(label, meta)
#     self._resolved = True
#     self._debug = debug
#     self._path = [self.meta()['path']]
#     self._resolved = False
#
#   def identify(self):
#     super().identify()
#     if self._debug:
#       print("Path:", self._path)
#
#   def datum_type(self):
#     return self.type()
#
#   def is_subtype_of(self, other):
#     return False
#
#
# class GenericFileFactory(InputFactory):
#   @classmethod
#   def type(cls) -> InputType:
#     return generic_file
#
#   @classmethod
#   def schema(cls):
#     return {
#       'schema': {
#         'path': {'type': 'string'},
#         'label': {'type': 'string'}
#       },
#       'nullable': True
#     }
#
#   @classmethod
#   def build(cls, label: str, meta: Dict, debug=False) -> GenericFile:
#     return GenericFile(label, meta)
