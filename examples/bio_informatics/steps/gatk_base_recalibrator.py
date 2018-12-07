from typing import List

from examples.bio_informatics.data_types.bed import Bed
from examples.bio_informatics.data_types.ref_fasta import RefFasta
from types.common_data_types import File, String, Array
from Tool.tool import Tool, ToolOutput, ToolInput


class GatkBaseRecalibrator(Tool):
    @staticmethod
    def tool():
        return "gatk-base-recalibrator"

    @staticmethod
    def base_command():
        return "javac"

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("inputBam_BaseRecalibrator", File()),
            ToolInput("outputfile_BaseRecalibrator", String()),
            ToolInput("reference", RefFasta()),
            ToolInput("known", Array(File())),
            ToolInput("bedFile", Bed()),
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("out", File())
        ]
