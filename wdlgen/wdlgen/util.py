from abc import ABC, abstractmethod


class WdlBase(ABC):
    @abstractmethod
    def get_string(self):
        raise Exception("Subclass must override .get_string() method")
