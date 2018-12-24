from abc import ABC

from Pipeline import CommandTool, ToolInput, ToolArgument


class Gatk3Base(CommandTool, ABC):

    @staticmethod
    def base_command():
        return ["java"]

    def inputs(self):
        return [

        ]

    @staticmethod
    @abstractmethod
    def docker():
        raise Exception("An error likely occurred when resolving the method order for docker for the Gatk3 classes "
                        "or you're trying to execute the docker method of the base class (ie, don't do that). "
                        "The method order resolution must preference Gatkbase subclasses, "
                        "and the subclass must contain a definition for docker.")

    def arguments(self):
        return [
            ToolArgument("./test/test-files", position=2, prefix="-Djava.io.tmpdir=",
                         separate_value_from_prefix=False),
            ToolArgument("/usr/GenomeAnalysisTK.jar", position=3, prefix="-jar"),
        ]