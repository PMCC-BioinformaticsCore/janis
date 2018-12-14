from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.dbsnp import Dbsnp
from Pipeline.bioinformatics.data_types.fasta import Fasta
from Pipeline.bioinformatics.data_types.vcfidx import VcfIdx
from Pipeline import String, Array, File, CommandTool, ToolOutput, ToolInput, ToolArgument


# class GatkMutect(CommandTool):
#     tumor = ToolInput("tumor", Bam())
#     normal = ToolInput("normal", Bam())
#     reference = ToolInput("reference", Fasta())
#     bedFile = ToolInput("bedFile", Bed())
#     outputFilename = ToolInput("outputFilename", String())
#     dbsnp = ToolInput("dbsnp", Dbsnp())
#     cosmic = ToolInput("cosmic", Array(File()))
#
#     vcf = ToolOutput("vcf", VcfIdx())
#
#     @staticmethod
#     def tool():
#         return "gatk-mutect2"
#
#     @staticmethod
#     def base_command():
#         return "javac"


class Gatkmutect2(CommandTool):
    inputBam_tumor = ToolInput("inputBam_tumor", File(), position=5, prefix="-I:Tumor")
    inputBam_normal = ToolInput("inputBam_normal", File(), position=6, prefix="-I:Normal")
    bedFile = ToolInput("bedFile", File(), position=7, prefix="-L")
    reference = ToolInput("reference", File(), position=8, prefix="-R")
    dbsnp = ToolInput("dbsnp", Array(File()), position=10, prefix="--dbsnp")
    cosmic = ToolInput("cosmic", Array(File()), position=11, prefix="--cosmic")

    output_mutect2 = ToolOutput("output_mutect2", File(), glob='$(inputs.outputfile_name)')

    @staticmethod
    def tool():
        return "GatkMutect2"

    @staticmethod
    def base_command():
        return ['java']

    @staticmethod
    def docker():
        # TODO: Investigate why this is different to other gatk broad dockers
        return "scidap/gatk:v3.5"

    @staticmethod
    def doc():
        return None

    def arguments(self):
        return [ToolArgument("/usr/GenomeAnalysisTK.jar", position=3, prefix="-jar"),
                ToolArgument("MuTect2", position=4, prefix="-T")]

    java_arg = ToolInput("java_arg", String(optional=True), position=1, default="-Xmx4g")
    outputfile_name = ToolInput("outputfile_name", String(optional=True), position=9, prefix="-o")


if __name__ == "__main__":
    print(Gatkmutect2().help())
