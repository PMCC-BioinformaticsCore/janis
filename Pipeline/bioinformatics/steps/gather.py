from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline import File, Tool, ToolOutput, ToolInput


class Gather(Tool):
    bamFile = ToolInput("bamFile", Bam())
    bamIndex = ToolInput("bamIndex", File())

    out = ToolOutput("out", Bam())

    @staticmethod
    def tool():
        return "gather"

    @staticmethod
    def base_command():
        return "javac"
