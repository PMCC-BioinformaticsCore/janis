import glob

from pipeline_definition.types.input_type import InputFactory
from pipeline_definition.types.input_type import Input


class BamFactory(InputFactory):
  @classmethod
  def type(cls):
    return 'bam'

  @classmethod
  def label(cls):
    return 'BAM file'

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
  def build(cls, input_dict, debug=False):
    return BamInput(input_dict, debug=debug)


class BamInput(Input):

  def __init__(self, input_dict, debug=False):
    super().__init__(input_dict)
    self.__debug = debug
    self.path = None
    self._files = []
    self._resolved = False

    if self.meta is not None:
      self.path = self.meta().get("path")

  def translate_for_workflow(self):
    raise Exception('translate_for_input not implemented for BAM files')

  def translate_for_input(self):
    if self._resolved:
      fd = [{'class': 'File', 'path': f} for f in self._files]
    else:
      fd = {'class': 'File', 'path': self.meta()['path']}

    return {self.id(): fd}

  def resolve(self):
    self._resolved = True
    self._files = glob.glob(self.meta()['path'])

  def identify(self):
    super().identify()
    if self.__debug:
      print("Path:", self.path)

  def datum_type(self):
    return self.type()

  def is_subtype_of(self, other):
    return False
