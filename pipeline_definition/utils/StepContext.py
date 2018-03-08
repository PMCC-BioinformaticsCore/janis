import json

class StepContext:
    def __init__(self, step):
        self.__step = step
        self.__branchOutputsStack = list()
        self.__dependecnyContexts = list()

    def inheritContextOfBranch(self, prevCtx):
        outputsStackFromPrevContext = prevCtx.outputStackOfBranch()
        if not outputsStackFromPrevContext:
            return
        self.__branchOutputsStack.extend(outputsStackFromPrevContext)

    def outputStackOfBranch(self):
        outputStack = list()

        providesDict = {}
        providesDict[self.__step.id()] = self.__step.provides()
        outputDict = {}
        outputDict[self.__step.tag()] = providesDict
        outputStack.append(outputDict)

        outputStack.extend(self.__branchOutputsStack)

        return outputStack

    def addDependencyContextFrom(self, stepCtx ):
        dependencyContext = stepCtx.provides()
        self.__dependecnyContexts.append(dependencyContext)

    def provides(self):
        providesDict = {}
        providesDict[self.__step.id()] = self.__step.provides()
        outputDict = {}
        outputDict[self.__step.tag()] = providesDict
        return outputDict

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__,sort_keys=True, indent=4)

    def print(self):
        print(self.toJSON())
