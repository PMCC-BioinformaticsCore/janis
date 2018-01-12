#
#
#
from abc import ABC, abstractmethod


class StepFactory(ABC):
    @classmethod
    @abstractmethod
    def build(cls, yml):
        pass

    @classmethod
    @abstractmethod
    def describe(cls):
        pass

    @classmethod
    @abstractmethod
    def type(cls):
        pass

    @classmethod
    @abstractmethod
    def description(cls):
        pass

    @classmethod
    @abstractmethod
    def label(cls):
        pass


class Step(ABC):
    @abstractmethod
    def __init__(self, yml):
        self.type = yml['type']

    def type(self):
        return self.type

    @abstractmethod
    def can_default(self):
        pass
