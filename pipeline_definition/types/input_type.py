#
# Input type base class
#
from abc import ABC, abstractmethod


class InputFactory(ABC):
    @classmethod
    @abstractmethod
    def type(cls):
        pass

    @classmethod
    @abstractmethod
    def label(cls):
        pass

    @classmethod
    @abstractmethod
    def description(cls):
        pass

    @classmethod
    @abstractmethod
    def describe(cls):
        pass

    @classmethod
    @abstractmethod
    def build(cls, meta):
        pass

    @classmethod
    def buildFrom(self, dict):
        type = self.type()
        print(type, "factory: Building from", dict)
        return self.build( dict )



class Input(ABC):
    def __init__(self, dict):
        self.label = next(iter(dict.keys()))

        meta = next(iter(dict.values()))
        self.type = next(iter(meta.keys()))
        self.meta = next(iter(meta.values()))

        self.identify()

    def identify(self):
        print("Instance: [", self.label, " - ", self.type, " - ", self.meta, " ]" )


    #@classmethod
    #@abstractmethod
    #def type(self):
    #    pass

    #@classmethod
    #@abstractmethod
    #def label(self):
    #    pass

    #@classmethod
    #@abstractmethod
    #def description(self):
    #    pass

    #@classmethod
    #@abstractmethod
    #def save(self):
    #    pass

    #@classmethod
    #@abstractmethod
    #def file_set(self):
    #    pass


# class InputSet(ABC):
#     @classmethod
#     @abstractmethod
#     def __init__(self, yml):
#         self.type = yml['type']
#
#     @classmethod
#     @abstractmethod
#     def input_set(self):
#         pass
#
#     @classmethod
#     def type(self):
#         return self.type
