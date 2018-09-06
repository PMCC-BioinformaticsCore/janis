#
# Input type base class
#
from abc import ABC, abstractmethod


class InputFactory(ABC):
  @classmethod
  @abstractmethod
  def type(cls):
    pass

  @classmethod
  @abstractmethod
  def label(cls):
    pass

  @classmethod
  @abstractmethod
  def description(cls):
    pass

  @classmethod
  @abstractmethod
  def describe(cls):
    pass

  @classmethod
  @abstractmethod
  def build(cls, input_dict):
    pass

  @classmethod
  def build_from(cls, input_dict):
    input_type = cls.type()
    print(input_type, "factory: Building from", input_dict)
    obj = cls.build(input_dict)
    obj.identify()
    return obj


class Input(ABC):
  def __init__(self, input_dict):
    self.__id = next(iter(input_dict.keys()))

    meta = next(iter(input_dict.values()))
    self.__type = next(iter(meta.keys()))
    self.__meta = next(iter(meta.values()))

  def identify(self):
    print("Instance: [", self.id, " - ", self.type, " - ", self.meta, " ]")

  def id(self):
    return self.__id

  def type(self):
    return self.__type

  def meta(self):
    return self.__meta

  @abstractmethod
  def datum_type(self):
    # A datum_type
    pass

  @abstractmethod
  def is_subtype_of(self, other):
    pass

  # @abstractmethod
  # def processMeta(self, meta):
  #    pass

  # @classmethod
  # @abstractmethod
  # def type(self):
  #    pass

  # @classmethod
  # @abstractmethod
  # def label(self):
  #    pass

  # @classmethod
  # @abstractmethod
  # def description(self):
  #    pass

  # @classmethod
  # @abstractmethod
  # def save(self):
  #    pass

  # @classmethod
  # @abstractmethod
  # def file_set(self):
  #    pass

# class InputSet(ABC):
#     @classmethod
#     @abstractmethod
#     def __init__(self, yml):
#         self.type = yml['type']
#
#     @classmethod
#     @abstractmethod
#     def input_set(self):
#         pass
#
#     @classmethod
#     def type(self):
#         return self.type
