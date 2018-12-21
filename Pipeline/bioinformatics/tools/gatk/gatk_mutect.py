from Pipeline import String, Array, File, CommandTool, ToolOutput, ToolInput, ToolArgument, Filename
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.vcfidx import VcfIdx


class GatkMutect2(CommandTool):
    inputBam_tumor = ToolInput("inputBam_tumor", BamPair(), position=5, prefix="-I:Tumor")
    inputBam_normal = ToolInput("inputBam_normal", BamPair(), position=6, prefix="-I:Normal")
    bedFile = ToolInput("bedFile", Bed(), position=7, prefix="-L")
    reference = ToolInput("reference", FastaWithDict(), position=8, prefix="-R")
    dbsnp = ToolInput("dbsnp", Array(VcfIdx()), position=10, prefix="--dbsnp")
    cosmic = ToolInput("cosmic", Array(VcfIdx()), position=11, prefix="--cosmic")

    output = ToolOutput("output", File(), glob='$(inputs.outputfile_name)')

    @staticmethod
    def tool():
        return "GatkMutect2"

    @staticmethod
    def base_command():
        return ['java']

    @staticmethod
    def docker():
        return "broadinstitute/gatk3:3.7-0"

    @staticmethod
    def doc():
        return None

    def arguments(self):
        return [
            ToolArgument("/usr/GenomeAnalysisTK.jar", position=3, prefix="-jar"), # ToolArgument("/usr/GenomeAnalysisTK.jar", position=3, prefix="-jar"),
            ToolArgument("MuTect2", position=4, prefix="-T")
        ]

    java_arg = ToolInput("java_arg", String(optional=True), position=1, default="-Xmx4g")
    outputfile_name = ToolInput("outputfile_name", Filename(), position=9, prefix="-o")


if __name__ == "__main__":
    print(GatkMutect2().help())
