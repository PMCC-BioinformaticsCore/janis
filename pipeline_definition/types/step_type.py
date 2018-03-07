#
#
#
from abc import ABC, abstractmethod

import sys

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
        self.__id = next(iter(dict.keys()))

        meta = next(iter(dict.values()))

        if meta is not None:
            self.__type = Step.selectTypeNameFrom(meta)
            self.__meta = meta[self.__type]
        else:
            self.__type = self.id
            self.__meta = None

        self.__tag = meta.get('tag')
        if self.__tag is None:
            self.__tag = "default"

        self.__gather = meta.get('gather')
        if self.__gather is not None:
            if isinstance(self.__gather, list):

                sortedList = sorted(self.__gather)
                ttag = None
                for item in sortedList:
                    if ttag is None:
                        ttag = item
                    else:
                        ttag = ttag + ":" + item

                self.__tag = ttag


    def tag(self):
        return self.__tag

    def gather(self):
        return self.__gather

    def id(self):
        return self.__id

    @abstractmethod
    def provides(self):
        # A set of optionally tagged input data
        pass

    @abstractmethod
    def requires(self):
        # A set of optionally tagged output data
        pass

    def identify(self):
        print("Instance: [", self.__id, " - ", self.__type, " - ", self.__meta, " ]" )

    @staticmethod
    def selectTypeNameFrom( meta ):
        selection = None
        for candidate in iter(meta.keys()):
            if candidate == 'tag':
                continue
            if candidate == 'gather':
                continue
            selection = candidate
            break
        return selection


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

class TaggedDatum(ABC):
    @abstractmethod
    def tags(self):
        # A set tags that can select among similar types
        pass

    @abstractmethod
    def datum_type(self):
        # A datum_type
        pass

    def satisfies(self, datum):
        # A concrete implementation here
        pass

