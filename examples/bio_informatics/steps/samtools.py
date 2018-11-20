from typing import List

from examples.bio_informatics.data_types.bam import Bam
from examples.bio_informatics.data_types.sam import Sam
from pipeline_definition.types.common_data_types import String
from pipeline_definition.types.tool import Tool, ToolOutput, ToolInput


class SamTools(Tool):
    @staticmethod
    def tool():
        return "samtools-view"

    @staticmethod
    def supported_translations() -> List[str]:
        return ["cwl"]

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("input", Sam()),
            ToolInput("outputName", String())
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("out", Bam())
        ]