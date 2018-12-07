from typing import List

from examples.bio_informatics.data_types.bam import Bam
from examples.bio_informatics.data_types.sam import Sam
from types.common_data_types import String
from Tool.tool import Tool, ToolOutput, ToolInput


class SamTools(Tool):
    @staticmethod
    def tool():
        return "samtools-view"

    @staticmethod
    def base_command():
        return "javac"

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("input", Sam()),
            ToolInput("outputName", String())
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("out", Bam())
        ]