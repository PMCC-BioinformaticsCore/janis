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
    CHROMOSOME = "chromosome"
    THIRTYX = "30x"
    NINETYX = "90x"
    THREEHUNDREDX = "300x"
    ALL = [TARGETED, EXOME, CHROMOSOME, THIRTYX, NINETYX, THREEHUNDREDX]

    @staticmethod
    def key(): return CaptureType.KEY

    @staticmethod
    def symbols(): return CaptureType.ALL


HINTS: List[Type[Hint]] = [CaptureType]
