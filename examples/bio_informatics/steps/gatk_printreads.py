from typing import List

from examples.bio_informatics.data_types.bam import Bam
from examples.bio_informatics.data_types.bed import Bed
from examples.bio_informatics.data_types.ref_fasta import RefFasta
from pipeline_definition.types.common_data_types import File, String
from pipeline_definition.types.tool import Tool, ToolOutput, ToolInput


class GatkPrintReads(Tool):
    @staticmethod
    def tool():
        return "gatk-printreads"

    @staticmethod
    def supported_translations() -> List[str]:
        pass

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("inputBam_printReads", Bam()),
            ToolInput("input_baseRecalibrator", File()),
            ToolInput("reference", RefFasta()),
            ToolInput("outputfile_printReads", String()),
            ToolInput("bedFile", Bed())
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("out", File()),
            ToolOutput("out-idx", File())
        ]