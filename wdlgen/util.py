from abc import ABC, abstractmethod


class WdlBase(ABC):
    @abstractmethod
    def wdl(self):
        raise Exception("Subclass must override .wdl() method")
