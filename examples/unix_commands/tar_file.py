from pipeline_definition.types.input_type import InputFactory
from pipeline_definition.types.input_type import Input

import glob


class TarFileFactory(InputFactory):
  @classmethod
  def type(cls):
    return 'tar'

  @classmethod
  def label(cls):
    return 'TAR file'

  @classmethod
  def description(cls):
    return 'Compressed or uncompressed tar files.'

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
  def build(cls, input_dict, debug=False):
    return TarFile(input_dict, debug=debug)


class TarFile(Input):
  def translate_for_input(self):
    if self._resolved:
      fd = [{'class': 'File', 'path': f} for f in self._files]
    else:
      fd = {'class': 'File', 'path': self.meta()['path']}

    return {self.id(): fd}

  def resolve(self):
    self._resolved = True
    self._files = glob.glob(self.meta()['path'])

  def __init__(self, input_dict, debug=False):
    super().__init__(input_dict)
    self.__debug = debug
    self.path = None
    self._files = []
    self._resolved = False

    if self.meta is not None:
      self.path = self.meta().get("path")

  def identify(self):
    super().identify()
    if self.__debug:
      print("Path:", self.path)

  def datum_type(self):
    return self.type()

  def is_subtype_of(self, other):
    return False
