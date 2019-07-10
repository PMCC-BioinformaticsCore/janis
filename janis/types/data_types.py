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
from typing import Any, List, Optional

import cwlgen as cwl
from wdlgen import WdlType

from janis.utils.logger import Logger

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
        raise Exception(
            f"Unhandled primitive type {t}, expected one of {', '.join(NativeTypes.all)}"
        )

    @staticmethod
    def map_to_wdl(t: NativeType):
        import wdlgen as wdl

        if t == NativeTypes.kBool:
            return wdl.PrimitiveType.kBoolean
        elif t == NativeTypes.kInt:
            return wdl.PrimitiveType.kInt

        elif (
            t == NativeTypes.kLong
            or t == NativeTypes.kFloat
            or t == NativeTypes.kDouble
        ):
            return wdl.PrimitiveType.kFloat
        elif t == NativeTypes.kStr:
            return wdl.PrimitiveType.kString
        elif t == NativeTypes.kFile:
            return wdl.PrimitiveType.kFile
        elif t == NativeTypes.kStdout:
            return wdl.PrimitiveType.kFile
        elif t == NativeTypes.kDirectory:
            Logger.log(
                "Using data_type 'Directory' for wdl, this requires cromwell>=37 and language=development"
            )
            return wdl.PrimitiveType.kDirectory
        elif t == NativeTypes.kArray:
            return wdl.ArrayType.kArray
        raise Exception(
            f"Unhandled primitive type {t}, expected one of {', '.join(NativeTypes.all)}"
        )

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
            return {"type": "Directory", "path": "path/to/file"}
        elif t == NativeTypes.kArray:
            return []


class DataType(ABC):
    def __init__(self, optional=False):
        self.optional = optional
        self.is_prim = NativeTypes.is_primitive(self.primitive())

    @staticmethod
    @abstractmethod
    def name() -> str:
        raise Exception("Subclass MUST override name field")

    @classmethod
    def __hash__(cls):
        return cls.name()

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

    # The following methods don't need to be overriden, but can be

    def id(self):
        if self.optional:
            return f"Optional<{self.name()}>"
        return self.name()

    @abstractmethod
    def validate_value(self, meta: Any, allow_null_if_not_optional: bool) -> bool:
        pass

    def identify(self):
        print(self.id())

    def can_receive_from(self, other, source_has_default=False) -> bool:
        """
        Can this class receive from $other, likely going to be type(a) == type(b)
        :param other:
        :param source_has_default: If the source has default, then we can return true even if the source is optional
        :return:
        """
        if not isinstance(other.received_type(), type(self.received_type())):
            return False
        if self.optional or source_has_default:
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

    def _question_mark_if_optional(self, has_default: bool = False):
        return "?" if self.optional or has_default else ""

    def cwl_type(self, has_default=False):
        tp = NativeTypes.map_to_cwl(self.primitive())
        return (
            [tp, "null"] if self.optional and not has_default else tp
        )  # and not has_default

    def map_cwl_type(self, parameter: cwl.Parameter) -> cwl.Parameter:
        if not NativeTypes.is_valid(self.primitive()):
            raise Exception(
                f"{self.id()} must declare its primitive as one of the NativeTypes "
                f"({', '.join(NativeTypes.all)})"
            )

        tp = NativeTypes.map_to_cwl(self.primitive())
        parameter.type = [tp, "null"] if self.optional else tp
        parameter.secondaryFiles = self.secondary_files()
        return parameter

    def cwl_input(self, value: Any):
        return value

    def wdl(self, has_default=False) -> WdlType:
        qm = self._question_mark_if_optional(has_default)
        return WdlType.parse_type(NativeTypes.map_to_wdl(self.primitive()) + qm)

    # def default(self):
    #     return self.default_value
