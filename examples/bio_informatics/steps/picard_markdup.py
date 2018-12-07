from typing import List

from types.common_data_types import Array, String, Number, File
from Tool.tool import Tool, ToolOutput, ToolInput


class PicardMarkDup(Tool):
    @staticmethod
    def tool():
        return "picard-markdups"

    @staticmethod
    def base_command():
        return "javac"

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("inputFileName_markDups", Array(String())),
            ToolInput("outputFileName_markDups", String()),
            ToolInput("metricsFile", String(optional=True)),
            ToolInput("picard_markdup_tmpdir", String(optional=True)),
            ToolInput("maxRecordsInRam", Number(optional=True)),
            ToolInput("validation_stringency", String(optional=True)),
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("out", File()),
            ToolOutput("out_idx", File()),
            ToolOutput("metrics", File())
        ]
