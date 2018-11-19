###################
# Implementations #
###################
from pipeline_definition.types.data_types import DataType
from pipeline_definition.types.type_registry import register_type


class String(DataType):

    @staticmethod
    def name():
        return "String"

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class Number(DataType):

    @staticmethod
    def name():
        return "Number"

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class Boolean(DataType):

    @staticmethod
    def name():
        return "Boolean"

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class Array(DataType):
    def __init__(self, t: DataType, optional=False):
        if not isinstance(t, DataType):
            raise Exception(f"Type t ({type(t)}) must be an instance of 'Type'")

        self.__t = t
        super().__init__(optional)

    @staticmethod
    def name():
        return f"Array"

    def can_receive_from(self, other):
        if not self.__t.can_receive_from(other):
            return False
        return super().can_receive_from(other)

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class File(DataType):

    @staticmethod
    def name():
        return "File"

    def input_field_from_input(self, meta):
        # WDL: "{workflowName}.label" = meta["path"}
        return {
            "class": "File",
            "path": meta["path"]
        }


register_type(String)
register_type(Number)
register_type(Boolean)
register_type(Array)
register_type(File)
