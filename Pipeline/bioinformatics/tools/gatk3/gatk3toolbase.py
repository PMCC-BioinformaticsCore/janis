from abc import ABC, abstractmethod

from Pipeline import CommandTool, ToolInput, ToolArgument, Boolean


class Gatk3ToolBase(CommandTool, ABC):

    @staticmethod
    def base_command():
        return ["java"]

    @staticmethod
    @abstractmethod
    def analysis_type():
        raise Exception("Tool must provide the '--analysis-type': A complete list of tools (sometimes also called "
                        "walkers because they 'walk' through the data to perform analyses) is available in the "
                        "online documentation.https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_engine_CommandLineGATK.php#--analysis_type")

    def inputs(self):
        return [
            ToolInput("showFullBamList", Boolean(optional=True), prefix="--showFullBamList",
                      doc="Emit list of input BAM/CRAM files to log"),
            # --use_jdk_deflater
            # --use_jdk_inflater
            # ToolInput("c")
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
            ToolArgument(self.analysis_type(), position=4, prefix="-T")
        ]