###################
# Implementations #
###################
from typing import Dict, Any
import wdlgen as wdl
import cwlgen as cwl

from janis.types.data_types import DataType, NativeTypes, NativeType


class String(DataType):

    @staticmethod
    def name():
        return "String"

    @staticmethod
    def primitive():
        return NativeTypes.kStr

    def doc(self):
        return "A string"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "string", "required": True}

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))

    def can_receive_from(self, other):
        if isinstance(other, Filename):
            return True
        return super().can_receive_from(other)


class Filename(String):

    def __init__(self, extension: str=None):
        """
        :param extension: with no '.' (dot)
        """
        self.extension = extension
        super().__init__(optional=True, default=self.generated_filename())

    @staticmethod
    def name() -> str:
        return "Filename"

    @staticmethod
    def primitive() -> NativeType:
        return NativeTypes.kStr

    def cwl_type(self):
        self.optional = False
        t = super().cwl_type()
        self.optional = True
        return t

    def doc(self) -> str:
        return """
This class is a placeholder for generated filenames, by default it is optional and CAN be overrided, 
however the program has been structured in a way such that these names will be generated based on the step label. 
These should only be used when the tool _requires_ a filename to output and you aren't 
concerned what the filename should be. The Filename DataType should NOT be used as an output.
""".strip()

    @classmethod
    def schema(cls) -> Dict:
        pass

    def map_cwl_type(self, parameter: cwl.Parameter):
        super().map_cwl_type(parameter)
        parameter.default = self.generated_filename()

    def generated_filename(self, prefix: str=None) -> str:
        import uuid
        pre = (prefix + "-") if prefix is not None else ""
        ex = "" if self.extension is None else self.extension
        return pre + "generated-" + str(uuid.uuid1()) + ex

    def can_receive_from(self, other: DataType):
        # Specific override because Filename should be able to receive from string
        if isinstance(other, String):
            return True  # Always provides default, and is always optional
        return super().can_receive_from(other)


class Int(DataType):

    @staticmethod
    def name():
        return "Integer"

    @staticmethod
    def primitive():
        return NativeTypes.kInt

    def doc(self):
        return "An integer"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "number", "required": True}

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class Float(DataType):

    @staticmethod
    def name():
        return "Float"

    @staticmethod
    def primitive():
        return NativeTypes.kFloat

    def doc(self):
        return "A float"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "number", "required": True}

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class Double(DataType):

    @staticmethod
    def name():
        return "Double"

    @staticmethod
    def primitive():
        return NativeTypes.kDouble

    def doc(self):
        return "An integer"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "number", "required": True}

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class Boolean(DataType):

    @staticmethod
    def name():
        return "Boolean"

    @staticmethod
    def primitive():
        return NativeTypes.kBool

    def doc(self):
        return "A boolean"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "boolean", "required": True}

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class File(DataType):

    @staticmethod
    def name():
        return "File"

    @staticmethod
    def primitive():
        return NativeTypes.kFile

    def doc(self):
        return "A local file"

    @classmethod
    def schema(cls) -> Dict:
        return {
            "path": {"type": "string", "required": True}
        }

    def get_value_from_meta(self, meta):
        return meta.get("path")

    @staticmethod
    def cwl_input(value: Any):
        return {
            "class": cwl.CwlTypes.FILE,
            "path": value
        }

    # def wdl(self):
    #     # Todo: SECONDARY FILES
    #     return wdlgen.WdlType.parse_type(NativeTypes.map_to_wdl(self.primitive()) + self._question_mark_if_optional())


class Directory(DataType):
    @staticmethod
    def name():
        return "Directory"

    @staticmethod
    def primitive():
        return NativeTypes.kDirectory

    def doc(self):
        return "A directory of files"

    def get_value_from_meta(self, meta):
        return meta["path"]

    @classmethod
    def schema(cls) -> Dict:
        return {
            "path": {"type": "string", "required": True}
        }

    def input_field_from_input(self, meta):
        return meta["path"]

    @staticmethod
    def cwl_input(value: Any):
        # WDL: "{workflowName}.label" = meta["path"}
        return {
            "class": cwl.CwlTypes.DIRECTORY,
            "path": value
        }


class Array(DataType):

    def __init__(self, t: DataType, optional=False):
        if not isinstance(t, DataType):
            raise Exception(f"Type t ({type(t)}) must be an instance of 'Type'")

        self._t = t
        super().__init__(optional)

    def subtype(self):
        return self._t

    @staticmethod
    def name():
        return "Array"

    @staticmethod
    def primitive():
        return NativeTypes.kArray

    def id(self):
        if self._t is None:
            return super().id()
        t = self._t
        typed = f"Array<{t.id()}>"
        if self.optional:
            return f"Optional<{typed}>"
        return typed

    def doc(self):
        return "An array"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "array"}

    def cwl_type(self):
        inp = cwl.CommandInputArraySchema(
            items=self._t.cwl_type(),
            # label=None,
            # input_binding=None
        )
        return [inp, "null"] if self.optional else inp

    def map_cwl_type(self, parameter: cwl.Parameter) -> cwl.Parameter:
        parameter.type = cwl.CommandInputArraySchema(
            items=None,
            label=None,
            input_binding=None
        )
        return parameter

    def wdl(self) -> wdl.ArrayType:
        return wdl.ArrayType(self._t.wdl(), requires_multiple=False)

    def can_receive_from(self, other):
        if isinstance(other, Array):
            return self._t.can_receive_from(other._t)
        if not self._t.can_receive_from(other):
            return False
        return super().can_receive_from(other)

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class Stdout(File):

    @staticmethod
    def name():
        return "Stdout"

    @staticmethod
    def primitive():
        return NativeTypes.kStdout

    def __init__(self, subtype=File(), stdoutname=None):
        super().__init__(optional=False)
        self.subtype = subtype if subtype is not None else File()
        self.stdoutname = stdoutname

    def received_type(self):
        return self.subtype
