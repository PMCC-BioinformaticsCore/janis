###################
# Implementations #
###################
from typing import Dict, Any

from pipeline_definition.types.data_types import DataType
from pipeline_definition.types.type_registry import register_type


class String(DataType):

    @staticmethod
    def name():
        return "String"

    @staticmethod
    def primitive():
        return "string"

    @staticmethod
    def doc():
        return "A string"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "string", "required": True}

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class Number(DataType):

    @staticmethod
    def name():
        return "Number"

    @staticmethod
    def primitive():
        return "int"

    @staticmethod
    def doc():
        return "A number"

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
        return "boolean"

    @staticmethod
    def doc():
        return "A number"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "boolean", "required": True}

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))

    def cwl(self) -> Dict[str, Any]:
        return {
            **super().cwl(),
            "type": "Boolean",
        }


class File(DataType):

    @staticmethod
    def name():
        return "File"

    @staticmethod
    def primitive():
        return "File"

    @staticmethod
    def doc():
        return "A local file"

    def get_value_from_meta(self, meta):
        return meta["path"]

    @classmethod
    def schema(cls) -> Dict:
        return {
            "path": {"type": "string", "required": True}
        }

    def input_field_from_input(self, meta):
        # WDL: "{workflowName}.label" = meta["path"}
        return {
            "class": "File",
            "path": meta["path"]
        }


class Array(DataType):

    def __init__(self, t: DataType=File(), optional=False):
        if not isinstance(t, DataType):
            raise Exception(f"Type t ({type(t)}) must be an instance of 'Type'")

        self.__t = t
        super().__init__(optional)

    @staticmethod
    def name():
        return f"Array"


    @staticmethod
    def primitive():
        return "array"

    def id(self):
        if self.__t is None:
            return super().id()
        t = self.__t
        typed = f"Array<{t.id()}>"
        if self.optional:
            return f"Optional<{typed}>"
        return typed

    @staticmethod
    def doc():
        return "An array"

    @classmethod
    def schema(cls) -> Dict:
        return {"type": "array"}

    def cwl(self) -> Dict[str, Any]:
        return {
            "type": "Array",
        }

    def can_receive_from(self, other):
        if not self.__t.can_receive_from(other):
            return False
        return super().can_receive_from(other)

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


register_type(String)
register_type(Number)
register_type(Boolean)
register_type(Array)
register_type(File)
