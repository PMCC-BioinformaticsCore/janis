# process for adding new hints!
from abc import ABC, abstractmethod
from typing import List, Type


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

    TARGETED = "targeted"
    EXOME = "exome"
    CHROMOSOME = "chromosome"
    THIRTYX = "30x"
    NINETYX = "90x"
    THREEHUNDREDX = "300x"

    @staticmethod
    def key():
        return "captureType"

    @staticmethod
    def symbols():
        return [
            CaptureType.TARGETED,
            CaptureType.EXOME,
            CaptureType.CHROMOSOME,
            CaptureType.THIRTYX,
            CaptureType.NINETYX,
            CaptureType.THREEHUNDREDX,
        ]


class Engine(HintEnum):
    CROMWELL = "cromwell"

    @staticmethod
    def key():
        return "engine"

    @staticmethod
    def symbols():
        return [Engine.CROMWELL]


HINTS: List[Type[Hint]] = [CaptureType, Engine]
