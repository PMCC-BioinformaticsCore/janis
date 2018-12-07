from typing import List

from examples.bio_informatics.data_types.bam import Bam
from types.common_data_types import File, String, Number
from Tool.tool import Tool, ToolOutput, ToolInput


class PicardSortSam(Tool):
    @staticmethod
    def tool():
        return "picard-sortsam"

    @staticmethod
    def base_command():
        return "javac"

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("inputFileName_sortSam", File()),
            ToolInput("outputFileName_sortSam", String()),
            ToolInput("tmpdir", String()),
            ToolInput("validation_stringency", String()),
            ToolInput("maxRecordsInRam", Number(optional=True))
        ]

    def outputs(self) -> List[ToolOutput]:
        # Todo: Somehow join the outputs into one BAM/BAI file
        return [
            ToolOutput("out", Bam()),       # Bam file
            ToolOutput("indexes", File())   # Bai Index
        ]
