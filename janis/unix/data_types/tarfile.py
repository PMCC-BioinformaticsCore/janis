from janis.types.common_data_types import File


class TarFile(File):
    @staticmethod
    def name():
        return "TarFile"

    def doc(self):
        return "A tarfile, ending with .tar"


# tar_file = InputType('tar_file', label='a unix tar archive', description='A unix compressed or uncompressed tar archive')
#
#
# class TarFileFactory(InputFactory):
#   @classmethod
#   def type(cls) -> InputType:
#     return tar_file
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
#   def build(cls, label: str, input_dict: dict, debug=False):
#     return TarFile(label, input_dict, debug=debug)
#
#
# class TarFile(Input):
#   def translate_for_workflow(self) -> Dict[str, str]:
#     return {self.id(): 'File'}
#
#   def translate_for_input(self):
#     fd = {'class': 'File', 'path': self._path}
#     return {self.id(): fd}
#
#   def resolve(self):
#     pass
#
#   def __init__(self, label: str, input_dict: dict, debug=False):
#     super().__init__(label, input_dict)
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

