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
  def __init__(self, input_dict, debug=False):
    super().__init__(input_dict, debug)
    self.forwardPattern = None
    self.backwardPattern = None
    self.__debug = debug

    if self.meta is not None:
      self.forwardPattern = self.meta().get("forward-pattern")
      self.backwardPattern = self.meta().get("backward-pattern")

  def identify(self):
    super().identify()
    if self.__debug:
      print("Forward Pattern:", self.forwardPattern)
      print("Backward Pattern:", self.backwardPattern)

  def datum_type(self):
    return self.type()

  def is_subtype_of(self, other):
    return False
