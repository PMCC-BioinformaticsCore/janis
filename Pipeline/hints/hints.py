# process for adding new hints!
from abc import ABC, abstractmethod
from typing import Union, List, Type


class Hint(ABC):
    # string inputs
    @staticmethod
    @abstractmethod
    def key():
        raise Exception("Must implement key() method")


class HintEnum(Hint, ABC):
    @staticmethod
    @abstractmethod
    def symbols():
        raise Exception("Must implement symbols() method")


class HintArray(Hint, ABC):
    @staticmethod
    @abstractmethod
    def items():
        raise Exception("Must implement items() method")


class CaptureType(HintEnum):

    KEY = "captureType"
    TARGETED = "targeted"
    EXOME = "exome"
    THIRTYX = "30x"
    NINETYX = "90x"
    THREEHUNDREDX = "300x"
    ALL = [TARGETED, EXOME, THIRTYX, NINETYX, THREEHUNDREDX]

    @staticmethod
    def key(): return CaptureType.KEY
    @staticmethod
    def symbols(): return CaptureType.ALL

HINTS: List[Type[Hint]] = [CaptureType]

def get_cwl_schema_for_recognised_hints():
    import cwlgen.cwlgen as cwl

    schema = cwl.CommandInputRecordSchema("hints")

    def prepare_hint(hint_class: Type[Hint]):
        name = hint_class.key() + "_schema"

        if issubclass(hint_class, HintEnum):
            # assume is an enum
            return cwl.CommandInputEnumSchema(label=hint_class.key(), name=name, symbols=hint_class.symbols())
        elif issubclass(hint_class, HintArray):
            return cwl.CommandInputArraySchema(items=hint_class.items(), label=hint_class.key())
        else:
            return "string?"

    schema.fields = [
        cwl.CommandInputRecordSchema.CommandInputRecordField(hint.key(), ["null", prepare_hint(hint)])
        for hint in HINTS]

    return schema
