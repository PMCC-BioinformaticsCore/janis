import json
from pipeline_definition.types.step_type import Step

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

    def providesFor(self, step):
        provides = self.provides()

        if not provides:
            return []

        tagSpecific = provides.get(step.tag())

        if not tagSpecific:
            return []

        stepSpecific = tagSpecific.get(step.id())
        if not stepSpecific:
            return []

        return stepSpecific

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__,sort_keys=True, indent=4)

    def print(self):
        print(self.toJSON())


    def mapInput(self, input):

        doc = {}
        #mapping provided?
        providedMapping = self.__step.providedValueForRequirement( input[Step.STR_ID] )
        if providedMapping:
            doc['provided'] = providedMapping
        else:
            doc['provided'] = ""

        candidates = {}

        if providedMapping:
            dependencySpec = Step.dependencySpecFrom(providedMapping)
            candidates = self.findMatchForInputDependency( input, dependencySpec )
        else:
            candidates = self.findMatchForInput(input)

        if not candidates:
            candidates['ERROR'] = "Failed to find any candidate!!!!!"

        doc['candidates'] = candidates
        return doc


    def findMatchForInputDependency(self, input, dependencySpec ):

        return self.findMatchForInput( input )

    def findMatchForInput(self, input):

        matches = {}
        pref = 1

        inputID = input[Step.STR_ID]
        inputType = input[Step.STR_TYPE]
        stepTag = self.__step.tag()

        for priorityEntry in self.__branchOutputsStack:
            matched = False
            tag, step = next(iter(priorityEntry.items()))
            stepName, outputs = next(iter(step.items()))

            for o in outputs:
                oID = o[Step.STR_ID]
                oType = o[Step.STR_TYPE]

                #Name and Type match is heighest priority - conclusive
                name = inputID
                if name == oID and inputType == oType:
                    matches[pref] = o
                    pref = pref + 1
                    break

                #If name starts with or equals to tag and type matches then take that
                name = stepTag
                if oID.startswith(name) and inputType == oType:
                    matches[pref] = o
                    pref = pref + 1
                    break

                if inputID == oID:
                    matches[pref] = o
                    pref = pref + 1
                    continue

                if inputType == oType:
                    matches[pref] = o
                    pref = pref + 1





        return matches