import json

class StepContext:
    def __init__(self, step):
        self.__step = step
        self.__availableOutputsStack = list()

    def inheritContextOf(self, prevCtx):
        outputsStackFromPrevContext = prevCtx.availableOutputStack()
        if not outputsStackFromPrevContext:
            return
        self.__availableOutputsStack.extend(outputsStackFromPrevContext)

    def availableOutputStack(self):
        outputStack = list()

        providesDict = {}
        providesDict[self.__step.id()] = self.__step.provides()
        outputDict = {}
        outputDict[self.__step.tag()] = providesDict
        outputStack.append(outputDict)

        outputStack.extend(self.__availableOutputsStack)

        return outputStack




    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__,sort_keys=True, indent=4)

    def print(self):
        print(self.toJSON())
