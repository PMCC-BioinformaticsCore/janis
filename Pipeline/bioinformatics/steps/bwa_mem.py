from Pipeline.bioinformatics.data_types.fastq import FastQ
from Pipeline.bioinformatics.data_types.sam import Sam
from Pipeline import String, Int, Array, Tool, ToolOutput, ToolInput


class BwaMem(Tool):
    reads = ToolInput("reads", Array())
    reference = ToolInput("reference", FastQ())
    outputFilename = ToolInput("outputFilename", String())
    readGroup = ToolInput("readGroup", String()),
    threads = ToolInput("threads", Int(optional=True)),
    min_std_max_min = ToolInput("min_std_max_min", Array(optional=True))

    out = ToolOutput("out", Sam())

    @staticmethod
    def tool():
        return "bwa-mem"

    @staticmethod
    def base_command():
        return "javac"

