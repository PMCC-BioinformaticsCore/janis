#
# Input type base class
#
from abc import ABC, abstractmethod

# The InputType defines a data object that is required as an input to a bit of software
# or produced by a bit of software. Typically, these will be an object in a store like a POSIX
# file system or they may be a collection of files or objects that go together. Examples are
# paired reads archives or indexed references.
#
# All data are output types, even if they are the output of an input declaration.
#
# This class represents the type. Steps can require and provide types. These work with
# - InputFactory classes which read input descriptions and provide
# - Input objects which represent the actual object in the store

# Type name quoted here because of Python's inability to handle circular dependencies
# https://www.python.org/dev/peps/pep-0484/#forward-references
from pipeline_definition.utils.errors import NotFoundException
from pipeline_definition.utils.registry import Registry

__input_types = Registry['InputType']()


def register_input_type(input_type: 'InputType'):
    __input_types.register(input_type.type_name(), input_type)


def get_input_type(type_name: str) -> 'InputType':
    try:
        return __input_types.get(type_name)
    except KeyError:
        raise NotFoundException(f'Input type {type_name} is not recognized. ' +
                                'This might mean a typo in the pipeline for a missing import.')


class InputType:
    def __init__(self, type_name: str, label: str = None, description: str = None):
        self._type_name = type_name
        self._label = label
        if description is None:
            self._description = label
        else:
            self._description = description
        register_input_type(self)

    def type_name(self) -> str:
        return self._type_name

    def label(self) -> str:
        return self._label

    def description(self) -> str:
        return self._description


class Input(ABC):
    def __init__(self, label: str, meta: dict):

        self.__id: str = label
        self.__type: str = meta["type"]
        self.__meta: dict = meta
        self.__debug = False

    @staticmethod
    def _get_type(type_name: str) -> InputType:
        return get_input_type(type_name)

    def identify(self):
        if self.__debug:
            print(f"Instance: [{self.id} - {self.type} - {self.meta}]")

    def id(self) -> str:
        # The id by which this input will be referred.
        return self.__id

    def type(self) -> InputType:
        # A string identifying the file type
        return self.__type

    def meta(self) -> dict:
        # Internal metadata required by this object
        return self.__meta

    @abstractmethod
    def datum_type(self):
        pass

    @abstractmethod
    def is_subtype_of(self, other) -> bool:
        pass

    def resolve(self):
        # Resolve actual object names in the appropriate store. For example, if the input is a query (regex)
        # the query is executed and resolved objects are returned as part of the translation. If the file
        # is a reference, its existence is checked.
        raise NotImplementedError("resolve not implemented for this type")

    def translate_for_input(self) -> dict:
        # - Generate the input object list. In the case of CWL, this is typically a separate yml file
        #   specifying a list of actual files.
        # - The expected return is a target language specific dictionary that the translator will render
        #   to the actual language text.
        raise NotImplementedError("A translation has been requested but has not been implemented for this input")

    def translate_for_workflow(self) -> dict:
        # - Generate the input stanza required by the target language. In the case of CWL, this is the 'inputs'
        #   section
        # - The expected return is a target language specific dictionary that the translator will render
        #   to the actual language text.
        raise NotImplementedError("A translation has been requested but has not been implemented for this input")


class InputFactory(ABC):
    @classmethod
    @abstractmethod
    def type(cls) -> InputType:
        # Input type that this factory generates.
        pass

    @classmethod
    def label(cls) -> str:
        # A brief description that can be rendered in UIs
        return cls.type().label()

    @classmethod
    def description(cls) -> str:
        # A description that can be rendered in UIs
        return cls.type().description()

    @classmethod
    @abstractmethod
    def schema(cls) -> dict:
        # The schema for this type.
        pass

    @classmethod
    @abstractmethod
    def build(cls, label: str, input_dict: dict) -> Input:
        # Build an Input object given the definition in the input_dict
        pass

    @classmethod
    def build_from(cls, label: str, input_dict: dict, debug: bool = False) -> Input:
        input_type = cls.type()
        if debug:
            print(input_type, "factory: Building from", input_dict)
        obj = cls.build(label, input_dict)
        obj.identify()
        return obj
