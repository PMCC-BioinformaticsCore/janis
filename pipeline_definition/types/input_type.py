#
# Input type base class
#
from abc import ABC, abstractmethod


class InputFactory(ABC):
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

    @classmethod
    @abstractmethod
    def emit(cls):
        pass


class InputSet(ABC):
    @abstractmethod
    def __init__(self, yml):
        self.type = yml['type']

    @abstractmethod
    def input_set(self):
        pass

    def type(self):
        return self.type


class Input(ABC):
    @abstractmethod
    def file_set(self):
        pass

    @abstractmethod
    def type(self):
        pass

    @abstractmethod
    def save(self):
        pass

    @abstractmethod
    def label(self):
        pass

    @abstractmethod
    def description(self):
        pass