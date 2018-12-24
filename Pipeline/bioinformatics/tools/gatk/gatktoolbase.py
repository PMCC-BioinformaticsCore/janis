from abc import ABC, abstractmethod

from Pipeline import CommandTool


class GatkToolBase(CommandTool, ABC):
    @staticmethod
    @abstractmethod
    def tool():
        raise Exception("Must override the 'tool' method")

    @abstractmethod
    def inputs(self):
        raise Exception("Must override the 'inputs' method")

    @staticmethod
    def doc():
        raise Exception("Must override the 'doc' method")

    @abstractmethod
    def arguments(self):
        return []
