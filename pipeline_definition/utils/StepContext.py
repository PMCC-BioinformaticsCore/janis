import json

class StepContext:
    def __init__(self, step):
        self.__step = step
        self.__availableOutputsStack = list()

        providesDict = {}
        providesDict[step.id()] = step.provides()
        outputDict = {}
        outputDict[step.tag()] = providesDict
        self.__availableOutputsStack.append(outputDict)

    def inheritContextOf(self, prevCtx):
        outputsStackFromPrevContext = prevCtx.availableOutputStack()
        if not outputsStackFromPrevContext:
            return












        pass

    def availableOutputStack(self):
        return self.__availableOutputsStack




    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__,sort_keys=True, indent=4)

    def print(self):
        print(self.toJSON())
