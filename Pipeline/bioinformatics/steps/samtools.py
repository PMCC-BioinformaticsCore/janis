from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.sam import Sam
from Pipeline import String, Tool, ToolOutput, ToolInput


class SamTools(Tool):
    inp = ToolInput("input", Sam()),
    outputName = ToolInput("outputName", String())

    out = ToolOutput("out", Bam())

    @staticmethod
    def tool():
        return "samtools-view"

    @staticmethod
    def base_command():
        return "javac"
