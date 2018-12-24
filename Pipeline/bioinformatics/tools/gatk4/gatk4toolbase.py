from abc import ABC, abstractmethod

from Pipeline import CommandTool, ToolInput, Boolean, ToolArgument


class Gatk4ToolBase(CommandTool, ABC):

    @staticmethod
    def base_command():
        return ["java"]

    def inputs(self):
        return [
            ToolInput("pg-tag", Boolean(optional=True), prefix="--add-output-sam-program-record",
                      doc="If true, adds a PG tag to created SAM/BAM/CRAM files.")
        ]

    @staticmethod
    @abstractmethod
    def docker():
        raise Exception("An error likely occurred when resolving the method order for docker for the Gatk classes "
                        "or you're trying to execute the docker method of the base class (ie, don't do that). "
                        "The method order resolution must preference Gatkbase subclasses, "
                        "and the subclass must contain a definition for docker.")

    def arguments(self):
        return [
            ToolArgument("./test/test-files", position=2, prefix="-Djava.io.tmpdir=",
                         separate_value_from_prefix=False),
            ToolArgument("/usr/GenomeAnalysisTK.jar", position=3, prefix="-jar"),
        ]