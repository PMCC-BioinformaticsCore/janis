#
#
#
from abc import ABC, abstractmethod


class StepFactory(ABC):
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
        obj = self.build(dict)
        obj.identify()
        return obj


class Step(ABC):
    def __init__(self, dict):
        self.id = next(iter(dict.keys()))

        meta = next(iter(dict.values()))

        if meta is not None:
            self.type = next(iter(meta.keys()))
            self.meta = next(iter(meta.values()))
        else:
            self.type = self.id
            self.meta = None


    def identify(self):
        print("Instance: [", self.id, " - ", self.type, " - ", self.meta, " ]" )

    #@abstractmethod
    #def type(self):
    #    pass

    # @classmethod
    #@abstractmethod
    #def label(self):
    #    pass

    # @classmethod
    #@abstractmethod
    #def description(self):
    #    pass

    # @classmethod
    #@abstractmethod
    #def save(self):
    #    pass

    # @classmethod
    #@abstractmethod
    #def can_default(self):
    #    pass
