from abc import ABC, abstractmethod
from typing import List, Dict

from pipeline_definition.types.input_type import InputType, Input


class Step(ABC):
  STR_ID = "id"
  STR_TYPE = "type"

  def __init__(self, input_dict, debug=False):
    self.__id = next(iter(input_dict.keys()))
    self.__debug = debug

    step_meta = next(iter(input_dict.values()))

    if step_meta is not None:
      self.__type = Step.select_type_name_from(step_meta)
      self.__meta = step_meta[self.__type]
    else:
      self.__type = self.id
      self.__meta = None

    self.__tag = None
    if step_meta is not None:
      self.__tag = step_meta.get('tag')

    if self.__tag is None:
      self.__tag = Step.default_step_tag_name()

  def tag(self) -> str:
    return self.__tag

  def type(self) -> InputType:
    return self.__type

  def id(self) -> str:
    return self.__id

  def meta(self) -> dict:
    return self.__meta

  def identify(self):
    if self.__debug:
      print("Instance: [", self.__id, " - ", self.__type, " - ", self.__meta, " ]")

  def provided_value_for_requirement(self, step_input: InputType):

    if self.__meta is None:
      return None

    provided = self.__meta.get(step_input.type_name())
    return provided

  @abstractmethod
  def provides(self) -> Dict[str, InputType]:
    # A map of output ids and output types. For example, an aligner might provide
    # {'aligned_file': bam_file}. Dependent steps can then access this as
    # align_step_id/aligned_file
    raise RuntimeError("Please provide implementation")

  @abstractmethod
  def requires(self) -> List[InputType]:
    # Return the input types required by this step
    raise RuntimeError("Please provide implementation")

  @abstractmethod
  def translate(self, mapped_inputs) -> Dict[str, str]:
    # Return a language specific dictionary that will be translated
    # to the output text.
    raise RuntimeError("Please provide implementation")

  @abstractmethod
  def cores(self) -> int:
    # Number of CPUS
    raise RuntimeError("Please provide implementation")

  @abstractmethod
  def ram(self) -> int:
    # Amount of ram
    raise RuntimeError("Please provide implementation")

  @staticmethod
  def select_type_name_from(meta) -> str:
    selection = None
    for candidate in iter(meta.keys()):
      if candidate == 'tag':
        continue
      if candidate == 'input_scope':
        continue
      selection = candidate
      break
    return selection

  def validate_input_output_spec(self):
    pass

    # output_specs = self.provides()
    # if output_specs:
    #   for ospec in output_specs:
    #     if not isinstance(ospec, InputType):
    #       raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")
    #
    #     if Step.STR_ID not in ospec:
    #       raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")
    #     name = ospec.get(Step.STR_ID)
    #     if not name:
    #       raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")
    #
    #     if Step.STR_TYPE not in ospec:
    #       raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")
    #     type = ospec.get(Step.STR_TYPE)
    #     if not type:
    #       raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")

  @staticmethod
  def default_step_tag_name():
    return 'untagged'

  @staticmethod
  def dependency_spec_from(requirement_value):

    tag = None
    step = None
    output = None

    parts = requirement_value.split(".")

    head = parts[0]
    if head.startswith("#"):
      tag = head[1:]

    return {
      'tag': tag,
      'step': step,
      'output': output
    }


class StepFactory(ABC):
  @classmethod
  @abstractmethod
  def type(cls) -> str:
    pass

  @classmethod
  @abstractmethod
  def label(cls) -> str:
    pass

  @classmethod
  def description(cls) -> str:
    return cls.label()

  @classmethod
  @abstractmethod
  def schema(cls) -> dict:
    pass

  @classmethod
  @abstractmethod
  def build(cls, meta, debug=False) -> Step:
    pass

  @classmethod
  def support_translations(cls) -> List[str]:
    return ['cwl']

  @classmethod
  def build_from(cls, step_dict, debug=False) -> Input:
    step_type = cls.type()
    if debug:
      print(step_type, "factory: Building from", step_dict)
    obj = cls.build(step_dict, debug=debug)
    obj.identify()
    return obj


class TaggedDatum(ABC):
  @abstractmethod
  def tags(self):
    # A set tags that can select among similar types
    pass

  @abstractmethod
  def datum_type(self):
    # A datum_type
    pass

  def satisfies(self, datum):
    # A concrete implementation here
    pass
