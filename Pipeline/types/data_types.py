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

from Pipeline.translations.cwl.cwl import Cwl
from Pipeline.translations.wdl.Wdl import Wdl

NativeType = str

class NativeTypes:
    kStr: NativeType = "str"
    kInt: NativeType = "int"
    kLong: NativeType = "long"
    kFloat: NativeType = "float"
    kBool: NativeType = "bool"
    kDouble: NativeType = "double"
    kFile: NativeType = "file"
    kDirectory: NativeType = "dir"
    kArray: NativeType = "array"

    _primitives: List[NativeType] = [kStr, kInt, kLong, kFile, kBool, kDouble]
    all: List[NativeType] = _primitives + [kFile, kDirectory, kArray]

    @staticmethod
    def is_primitive(t: NativeType) -> bool:
        return t in NativeTypes._primitives

    @staticmethod
    def is_valid(t: str):
        return t in NativeTypes.all

    @staticmethod
    def map_to_cwl(t: NativeType):
        if t == NativeTypes.kBool:
            return Cwl.PRIMITIVES.kBOOLEAN
        elif t == NativeTypes.kInt:
            return Cwl.PRIMITIVES.kINT
        elif t == NativeTypes.kLong:
            return Cwl.PRIMITIVES.kLONG
        elif t == NativeTypes.kFloat:
            return Cwl.PRIMITIVES.kFLOAT
        elif t == NativeTypes.kDouble:
            return Cwl.PRIMITIVES.kDOUBLE
        elif t == NativeTypes.kStr:
            return Cwl.PRIMITIVES.kSTRING
        elif t == NativeTypes.kFile:
            return Cwl.PRIMITIVES.kFILE
        elif t == NativeTypes.kDirectory:
            return Cwl.PRIMITIVES.kDIRECTORY
        elif t == NativeTypes.kArray:
            return Cwl.PRIMITIVES.kARRAY
        raise Exception(f"Unhandled primitive type {t}, expected one of {', '.join(NativeTypes.all)}")

    @staticmethod
    def map_to_wdl(t: NativeType):
        if t == NativeTypes.kBool:
            return Wdl.PRIMITIVES.kBOOLEAN
        elif t == NativeTypes.kInt:
            return Wdl.PRIMITIVES.kINT

        elif t == NativeTypes.kLong or t == NativeTypes.kFloat or t == NativeTypes.kDouble:
            return Wdl.PRIMITIVES.kFLOAT
        elif t == NativeTypes.kStr:
            return Wdl.PRIMITIVES.kSTRING
        elif t == NativeTypes.kFile:
            return Wdl.PRIMITIVES.kFILE
        elif t == NativeTypes.kDirectory:
            return Wdl.PRIMITIVES.kINT
        elif t == NativeTypes.kArray:
            return Wdl.PRIMITIVES.kARRAY
        raise Exception(f"Unhandled primitive type {t}, expected one of {', '.join(NativeTypes.all)}")

    @staticmethod
    def default_value(t: NativeType):
        if t == NativeTypes.kBool:
            return True
        elif t == NativeTypes.kInt:
            return 0
        elif t == NativeTypes.kLong:
            return 0.0
        elif t == NativeTypes.kFloat:
            return 0.0
        elif t == NativeTypes.kDouble:
            return 0.0
        elif t == NativeTypes.kStr:
            return "nothing"
        elif t == NativeTypes.kFile:
            return {"type": "File", "path": "path/to/file"}
        elif t == NativeTypes.kDirectory:
            return { "type": "Directory", "path": "path/to/file"}

class DataType(ABC):

    def __init__(self, optional=False):
        self.optional = optional
        self.is_prim = NativeTypes.is_primitive(self.primitive())

    @staticmethod
    @abstractmethod
    def name() -> str:
        raise Exception("Subclass MUST override name field")

    @staticmethod
    def secondary_files() -> Optional[List[str]]:
        return None

    @staticmethod
    @abstractmethod
    def primitive() -> NativeType:
        raise Exception("Subclass MUST override the 'schema' method")

    @staticmethod
    @abstractmethod
    def doc() -> str:
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
        if not NativeTypes.is_valid(self.primitive()):
            raise Exception(f"{self.id()} must declare its primitive as one of the NativeTypes "
                            f"({', '.join(NativeTypes.all)})")
        d = {
            "type": NativeTypes.map_to_cwl(self.primitive()) + self._question_mark_if_optional()
        }

        if self.doc():
            d["doc"] = self.doc()
        if self.secondary_files():
            d["secondaryFiles"] = self.secondary_files()

        return d

    def wdl(self):
        return NativeTypes.map_to_wdl(self.primitive())

