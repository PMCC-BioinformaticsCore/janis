
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.ref_fasta import RefFasta
from Pipeline import File, String, Tool, ToolOutput, ToolInput


class GatkPrintReads(Tool):
    inputBam_printReads = ToolInput("inputBam_printReads", Bam())
    input_baseRecalibrator = ToolInput("input_baseRecalibrator", File())
    reference = ToolInput("reference", RefFasta())
    outputfile_printReads = ToolInput("outputfile_printReads", String())
    bedFile = ToolInput("bedFile", Bed())

    out = ToolOutput("out", File())
    out_idx = ToolOutput("out_idx", File())

    @staticmethod
    def tool():
        return "gatk-printreads"

    @staticmethod
    def base_command():
        return "javac"