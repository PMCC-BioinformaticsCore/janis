from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.dbsnp import Dbsnp
from Pipeline.bioinformatics.data_types.fasta import Fasta
from Pipeline import String, Int, File, CommandTool, ToolOutput, ToolInput


class GatkHaplotypecaller(CommandTool):
    reference = ToolInput("reference", Fasta())
    outputfile_HaplotypeCaller = ToolInput("outputfile_HaplotypeCaller", String())
    dbsnp = ToolInput("dbsnp", Dbsnp())
    threads = ToolInput("threads", Int(optional=True))
    emitRefConfidence = ToolInput("emitRefConfidence", String())
    bedFile = ToolInput("bedFile", String())
    bamOutput = ToolInput("bamOutput", String(optional=True))

    out = ToolOutput("out", File())
    bamOut = ToolOutput("bamOut", Bam())

    @staticmethod
    def tool():
        return "gatk-haplotype"

    @staticmethod
    def base_command():
        return "javac"
