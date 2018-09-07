from pipeline_definition.types.input_type import InputFactory
from pipeline_definition.types.input_type import Input


class BAMFactory(InputFactory):
  @classmethod
  def type(cls):
    return 'BAM'

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
    return BAMInput(input_dict)


class BAMInput(Input):
  def __init__(self, input_dict, debug=False):
    super().__init__(input_dict)
    self.__debug = debug
    self.path = None

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
