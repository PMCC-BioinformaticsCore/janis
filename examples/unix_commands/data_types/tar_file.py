from pipeline_definition.types.input_type import InputFactory, InputType
from pipeline_definition.types.input_type import Input

tar_file = InputType('tar file', label='a unix tar archive', description='A unix compressed or uncompressed tar archive')


class TarFileFactory(InputFactory):
  @classmethod
  def type(cls) -> InputType:
    return tar_file

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
  def build(cls, input_dict, debug=False):
    return TarFile(input_dict, debug=debug)


class TarFile(Input):
  def translate_for_workflow(self) -> dict:
    raise Exception('Not yet implemented')

  def translate_for_input(self):
    fd = {'class': 'File', 'path': self._path}
    return {self.id(): fd}

  def resolve(self):
    pass

  def __init__(self, input_dict, debug=False):
    super().__init__(input_dict)
    self._resolved = True
    self._debug = debug
    self._path = [self.meta()['path']]
    self._resolved = False

  def identify(self):
    super().identify()
    if self._debug:
      print("Path:", self._path)

  def datum_type(self):
    return self.type()

  def is_subtype_of(self, other):
    return False
