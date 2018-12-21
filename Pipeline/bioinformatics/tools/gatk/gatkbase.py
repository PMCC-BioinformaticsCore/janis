from abc import ABC

from Pipeline import CommandTool, ToolInput, Boolean, ToolArgument


class GatkBase(CommandTool, ABC):

    @staticmethod
    def base_command():
        return ["java"]

    def inputs(self):
        return [
            ToolInput("pg-tag", Boolean(optional=True), prefix="--add-output-sam-program-record",
                      doc="If true, adds a PG tag to created SAM/BAM/CRAM files.")
        ]

    @staticmethod
    def docker():
        raise Exception("An error likely occurred when resolving the method order for docker for the Gatk classes. "
                        "The method order resolution must preference Gatkbase subclasses, "
                        "and the subclass must contain a definition for docker.")

    def arguments(self):
        return [
            ToolArgument("./test/test-files", position=2, prefix="-Djava.io.tmpdir=",
                         separate_value_from_prefix=False),
            ToolArgument("/usr/GenomeAnalysisTK.jar", position=3, prefix="-jar"),
        ]