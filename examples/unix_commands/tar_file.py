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
        'glob': {'type': 'string'},
        'label': {'type': 'string'}
      },
      'nullable': True
    }

  @classmethod
  def build(cls, input_dict, debug=False):
    return TarFile(input_dict, debug=debug)


class TarFile(Input):
  def resolve(self):
    self._files = glob.glob(self.__meta['glob'])

  def __init__(self, input_dict, debug=False):
    super().__init__(input_dict)
    self.__debug = debug
    self.path = None
    self._files = []

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
