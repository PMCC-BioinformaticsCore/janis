from typing import List

from examples.bio_informatics.data_types.bam import Bam
from examples.bio_informatics.data_types.bed import Bed
from examples.bio_informatics.data_types.dbsnp import Dbsnp
from examples.bio_informatics.data_types.ref_fasta import RefFasta
from examples.bio_informatics.data_types.vcfidx import VcfIdx
from types.common_data_types import String, Array, File
from Tool.tool import Tool, ToolOutput, ToolInput


class GatkMutect(Tool):
    @staticmethod
    def tool():
        return "gatk-mutect2"

    @staticmethod
    def base_command():
        return "javac"

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("tumor", Bam()),
            ToolInput("normal", Bam()),
            ToolInput("reference", RefFasta()),
            ToolInput("bedFile", Bed()),
            ToolInput("outputFilename", String()),
            ToolInput("dbsnp", Dbsnp()),
            ToolInput("cosmic", Array(File()))
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("vcf", VcfIdx())
        ]
