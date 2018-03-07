import json

class StepContext:
    def __init__(self, step):
        self.__step = step
        self.__prevStep = None
        self.__contextualInputSet = None
        self.__globalInputSet = None
        self.__globalOutputSet = None

    def setPrevStep( self, prevstep ):
        self.__prevStep = prevstep

    def setGlobalInputSet(self, globalInputSet):
        self.__globalInputSet = globalInputSet

    def setGlobalOutputSet(self, globalOutputSet):
        self.__globalOutputSet = globalOutputSet

    def setContextualInputSet(self, contextualInputSet):
        self.__contextualInputSet = contextualInputSet

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__,sort_keys=True, indent=4)

    def print(self):
        print(self.toJSON())
