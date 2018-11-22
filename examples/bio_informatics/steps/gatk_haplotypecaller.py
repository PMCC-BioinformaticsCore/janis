from typing import List

from examples.bio_informatics.data_types.bam import Bam
from examples.bio_informatics.data_types.dbsnp import Dbsnp
from examples.bio_informatics.data_types.ref_fasta import RefFasta
from pipeline_definition.types.common_data_types import String, Number, File
from pipeline_definition.types.tool import Tool, ToolOutput, ToolInput


class GatkHaplotypecaller(Tool):
    @staticmethod
    def tool():
        return "gatk-haplotype"

    @staticmethod
    def base_command():
        return "javac"

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("reference", RefFasta()),
            ToolInput("outputfile_HaplotypeCaller", String()),
            ToolInput("dbsnp", Dbsnp()),
            ToolInput("threads", Number(optional=True)),
            ToolInput("emitRefConfidence", String()),
            ToolInput("bedFile", String()),
            ToolInput("bamOutput", String(optional=True))
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("out", File()),
            ToolOutput("bamOut", Bam())
        ]
