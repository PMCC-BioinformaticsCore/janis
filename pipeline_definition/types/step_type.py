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
    STR_ID = "id"
    STR_TYPE = "type"

    def __init__(self, dict):
        self.__id = next(iter(dict.keys()))

        stepMeta = next(iter(dict.values()))

        if stepMeta is not None:
            self.__type = Step.selectTypeNameFrom(stepMeta)
            self.__meta = stepMeta[self.__type]
        else:
            self.__type = self.id
            self.__meta = None

        self.__tag = None
        if ( stepMeta is not None):
            self.__tag = stepMeta.get('tag')

        if self.__tag is None:
            self.__tag = Step.defaultSteptagName()

    def tag(self):
        return self.__tag

    def id(self):
        return self.__id

    def identify(self):
        print("Instance: [", self.__id, " - ", self.__type, " - ", self.__meta, " ]")

    def providedValueForRequirement(self, requirmentName):

        if self.__meta is None:
            return None

        provided = self.__meta.get(requirmentName)
        return provided

    @abstractmethod
    def provides(self):
        raise RuntimeError("Please provide implementation")

    @abstractmethod
    def requires(self):
        raise RuntimeError("Please provide implementation")

    @staticmethod
    def selectTypeNameFrom( meta ):
        selection = None
        for candidate in iter(meta.keys()):
            if candidate == 'tag':
                continue
            if candidate == 'input_scope':
                continue
            selection = candidate
            break
        return selection

    def validateInputOuputSpec(self):

        outputSpecs = self.provides()
        if outputSpecs:
            for ospec in outputSpecs:
                if not isinstance(ospec, dict ):
                    raise RuntimeError( "Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation." )


                if Step.STR_ID not in ospec:
                    raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")
                name = ospec.get(Step.STR_ID)
                if not name:
                    raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")

                if Step.STR_TYPE not in ospec:
                    raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")
                type = ospec.get(Step.STR_TYPE)
                if not type:
                    raise RuntimeError("Output spec provided by step " + self.id() + "[" + self.tag() + "] fails validation.")

    @staticmethod
    def defaultSteptagName():
        return 'untagged'

    @staticmethod
    def dependencySpecFrom( requirementValue ):

        tag = None
        step = None
        output = None

        parts = requirementValue.split(".")

        head = parts[0]
        if head.startswith("#"):
            tag = head[1:]

        return {
            'tag' : tag,
            'step' : step,
            'output' : output
        }

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

