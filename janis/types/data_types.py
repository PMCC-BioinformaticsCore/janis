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
from wdlgen import WdlType

import cwlgen as cwl

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
    kStdout: NativeType = "stdout"

    _primitives: List[NativeType] = [kStr, kInt, kFloat, kLong, kDouble, kBool, kDouble]
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
            return cwl.CwlTypes.BOOLEAN
        elif t == NativeTypes.kInt:
            return cwl.CwlTypes.INT
        elif t == NativeTypes.kLong:
            return cwl.CwlTypes.LONG
        elif t == NativeTypes.kFloat:
            return cwl.CwlTypes.FLOAT
        elif t == NativeTypes.kDouble:
            return cwl.CwlTypes.DOUBLE
        elif t == NativeTypes.kStr:
            return cwl.CwlTypes.STRING
        elif t == NativeTypes.kFile:
            return cwl.CwlTypes.FILE
        elif t == NativeTypes.kDirectory:
            return cwl.CwlTypes.DIRECTORY
        elif t == NativeTypes.kArray:
            return cwl.CwlTypes.ARRAY
        elif t == NativeTypes.kStdout:
            return cwl.CwlTypes.STDOUT
        raise Exception(f"Unhandled primitive type {t}, expected one of {', '.join(NativeTypes.all)}")

    @staticmethod
    def map_to_wdl(t: NativeType):
        import wdlgen as wdl
        if t == NativeTypes.kBool:
            return wdl.PrimitiveType.kBoolean
        elif t == NativeTypes.kInt:
            return wdl.PrimitiveType.kInt

        elif t == NativeTypes.kLong or t == NativeTypes.kFloat or t == NativeTypes.kDouble:
            return wdl.PrimitiveType.kFloat
        elif t == NativeTypes.kStr:
            return wdl.PrimitiveType.kString
        elif t == NativeTypes.kFile:
            return wdl.PrimitiveType.kFile
        elif t == NativeTypes.kStdout:
            return wdl.PrimitiveType.kFile
        elif t == NativeTypes.kDirectory:
            return None # wdl.PrimitiveType
        elif t == NativeTypes.kArray:
            return wdl.ArrayType.kArray
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
        elif t == NativeTypes.kArray:
            return []


class DataType(ABC):

    def __init__(self, optional=False, default=None):
        self.optional = optional
        self.is_prim = NativeTypes.is_primitive(self.primitive())
        self.default_value = default

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
        raise Exception("Subclass MUST override the 'primitive' method")

    @staticmethod
    @abstractmethod
    def doc() -> str:
        """
        Subclasses should override this class to provide additional information on how to
        correctly provide data to the class, what inputs it may have and what other types
        are compatible
        """
        raise Exception("Subclass MUST override the 'doc' field")

    @staticmethod
    def validate(meta: Any) -> bool:
        return True

    @classmethod
    # @abstractmethod
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
        if not isinstance(other.received_type(), type(self.received_type())):
            return False
        if self.optional:
            # If I'm optional I can receive from optional / non optional
            return True
        # If I'm not optional, I must receive from not optional
        return not other.optional

    def received_type(self):
        """
        The type that will be received if joined from this type
        :return: mostly self, except for STDOUT | STDERR | STDIN
        """
        return self

    def input_field_from_input(self, meta):
        """
        Method to convert the field definition into a generic CWL-esque response
        :param meta:
        :return:
        """
        return None

    def _question_mark_if_optional(self):
        return "?" if self.optional else ""

    def cwl_type(self):
        return NativeTypes.map_to_cwl(self.primitive()) + self._question_mark_if_optional()

    def map_cwl_type(self, parameter: cwl.Parameter) -> cwl.Parameter:
        if not NativeTypes.is_valid(self.primitive()):
            raise Exception(f"{self.id()} must declare its primitive as one of the NativeTypes "
                            f"({', '.join(NativeTypes.all)})")

        parameter.type = NativeTypes.map_to_cwl(self.primitive()) + self._question_mark_if_optional()
        parameter.secondaryFiles = self.secondary_files()
        return parameter

    @staticmethod
    def cwl_input(value: Any):
        return value

    def wdl(self) -> WdlType:
        return WdlType.parse_type(NativeTypes.map_to_wdl(self.primitive()) + self._question_mark_if_optional())

    def default(self):
        return self.default_value

