"""
    Custom Types::

    We want to define some standard interface that allows for the implementation of types.
    Types must be supported as Inputs and Outputs, and must be comparable, ie we can compare
    the output of one step to another step. We'll guarantee that

    We'll also may require separate interfaces for CWL and WDL implementations, we may not be
    able to genericise the interface enough to support an automatic conversion

    We are allowed to require that a type must register itself when loaded.

"""
from abc import ABC, abstractmethod
from typing import Any, List, Dict, Optional

allowed_primitives = ["null", "boolean", "int", "long", "float", "double", "string"]
allowed_types = allowed_primitives +  ["File", "Directory"]


class DataType(ABC):

    def __init__(self, optional=False):
        self.optional = optional
        self.is_prim = self.primitive() in allowed_primitives

    @staticmethod
    @abstractmethod
    def name():
        raise Exception("Subclass MUST override name field")

    @staticmethod
    def secondary_files() -> Optional[List[str]]:
        return None

    @staticmethod
    @abstractmethod
    def primitive():
        raise Exception("Subclass MUST override the 'schema' method")

    @staticmethod
    @abstractmethod
    def doc():
        """
        Subclasses should override this class to provide additional information on how to
        correctly provide data to the class, what inputs it may have and what other types
        are compatible
        """
        return None

    @staticmethod
    def validate(meta: Any) -> bool:
        return True

    @classmethod
    @abstractmethod
    def schema(cls) -> Dict:
        raise Exception("Subclass MUST override the 'schema' method")

    def get_value_from_meta(self, meta):
        return meta

    # The following methods don't need to be overriden, but can be

    def id(self):
        if self.optional:
            return f"Optional<{self.name()}>"
        return self.name()

    def identify(self):
        print(self.id())

    def can_receive_from(self, other) -> bool:
        """
        Can this class receive from $other, likely going to be type(a) == type(b)
        :param other:
        :return:
        """
        if not isinstance(other, type(self)):
            return False
        if self.optional:
            # If I'm optional I can receive from optional / non optional
            return True
        # If I'm not optional, I must receive from not optional
        return not other.optional

    def input_field_from_input(self, meta):
        """
        Method to convert the field definition into a generic CWL-esque response
        :param meta:
        :return:
        """
        return None

    def _question_mark_if_optional(self):
        return "?" if self.optional else ""

    def cwl(self) -> Dict[str, Any]:
        if self.primitive() not in allowed_types:
            raise Exception(f"{self.id()} must declare its primitive as one of {', '.join(allowed_types)}")
        d = {
            "type": self.primitive() + self._question_mark_if_optional()
        }

        if self.doc():
            d["doc"] = self.doc()
        if self.secondary_files():
            d["secondaryFiles"] = self.secondary_files()

        return d
