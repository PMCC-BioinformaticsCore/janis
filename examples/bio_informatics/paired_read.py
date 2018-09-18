import glob

from pipeline_definition.pipeline_translator import PipelineTranslatorException
from pipeline_definition.types.input_type import InputFactory
from pipeline_definition.types.input_type import Input


class PairedReadFactory(InputFactory):
  @classmethod
  def type(cls):
    return 'SequenceReadArchivePaired'

  @classmethod
  def label(cls):
    return 'paired read files'

  @classmethod
  def description(cls):
    return cls.label()

  @classmethod
  def describe(cls):
    return {
      'schema': {
        'label': {'type': 'string'},
        'forward-pattern': {'type': 'string', 'required': True},
        'backward-pattern': {'type': 'string'}
      },
      'nullable': True
    }

  @classmethod
  def build(cls, input_dict, debug=False):
    return PairedReadInput(input_dict, debug)


class PairedReadInput(Input):
  def translate(self):
    if self._resolved:
      forward = [{'class': 'File', 'path': f} for f in self.forward_files]
      backward = [{'class': 'File', 'path': f} for f in self.backward_files]
    else:
      forward = {'class': 'File', 'path': self.meta()['forward-pattern']}
      backward = {'class': 'File', 'path': self.meta()['backward-pattern']}

    return {
      self.id() + '_forward': forward,
      self.id() + '_backward': backward
    }

  def resolve(self):
    self.forward_files = glob.glob(self.forward_pattern)
    self.backward_files = glob.glob(self.backward_pattern)

    if len(self.forward_files) != len(self.backward_files):
      raise PipelineTranslatorException('The number of forward reads is different to the number of backward reads.')

    self._resolved = True

  def __init__(self, input_dict, debug=False):
    super().__init__(input_dict, debug)
    self._resolved = False
    self.forward_pattern = None
    self.backward_pattern = None
    self.__debug = debug
    self.forward_files = []
    self.backward_files = []

    if self.meta is not None:
      self.forward_pattern = self.meta().get("forward-pattern")
      self.backward_pattern = self.meta().get("backward-pattern")

  def identify(self):
    super().identify()
    if self.__debug:
      print("Forward Pattern:", self.forward_pattern)
      print("Backward Pattern:", self.backward_pattern)

  def datum_type(self):
    return self.type()

  def is_subtype_of(self, other):
    return False
