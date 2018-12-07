from typing import List

from examples.bio_informatics.data_types.bam import Bam
from types.common_data_types import File
from Tool.tool import Tool, ToolOutput, ToolInput


class Gather(Tool):
    @staticmethod
    def tool():
        return "gather"

    @staticmethod
    def base_command():
        return "javac"

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("bamFile", Bam()),
            ToolInput("bamIndex", File())
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("out", Bam())
        ]
