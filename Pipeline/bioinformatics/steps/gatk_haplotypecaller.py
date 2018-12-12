from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.dbsnp import Dbsnp
from Pipeline.bioinformatics.data_types.ref_fasta import RefFasta
from Pipeline import String, Int, File, Tool, ToolOutput, ToolInput


class GatkHaplotypecaller(Tool):
    reference = ToolInput("reference", RefFasta())
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
