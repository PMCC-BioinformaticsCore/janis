from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.dbsnp import Dbsnp
from Pipeline.bioinformatics.data_types.fasta import Fasta
from Pipeline.bioinformatics.data_types.vcfidx import VcfIdx
from Pipeline import String, Array, File, CommandTool, ToolOutput, ToolInput


class GatkMutect(CommandTool):
    tumor = ToolInput("tumor", Bam())
    normal = ToolInput("normal", Bam())
    reference = ToolInput("reference", Fasta())
    bedFile = ToolInput("bedFile", Bed())
    outputFilename = ToolInput("outputFilename", String())
    dbsnp = ToolInput("dbsnp", Dbsnp())
    cosmic = ToolInput("cosmic", Array(File()))

    vcf = ToolOutput("vcf", VcfIdx())

    @staticmethod
    def tool():
        return "gatk-mutect2"

    @staticmethod
    def base_command():
        return "javac"
