from typing import List

from Pipeline import Array, String, Int, File, Tool, ToolOutput, ToolInput


class PicardMarkDup(Tool):
    inputFileName_markDups = ToolInput("inputFileName_markDups", Array(String()))
    outputFileName_markDups = ToolInput("outputFileName_markDups", String())
    metricsFile = ToolInput("metricsFile", String(optional=True))
    picard_markdup_tmpdir = ToolInput("picard_markdup_tmpdir", String(optional=True))
    maxRecordsInRam = ToolInput("maxRecordsInRam", Int(optional=True))
    validation_stringency = ToolInput("validation_stringency", String(optional=True)),

    out = ToolOutput("out", File())
    out_idx = ToolOutput("out_idx", File())
    metrics = ToolOutput("metrics", File())


    @staticmethod
    def tool():
        return "picard-markdups"

    @staticmethod
    def base_command():
        return "javac"
