from typing import List

from examples.bio_informatics.data_types.bam import Bam
from examples.bio_informatics.data_types.bed import Bed
from examples.bio_informatics.data_types.dbsnp import Dbsnp
from examples.bio_informatics.data_types.ref_fasta import RefFasta
from examples.bio_informatics.data_types.vcfidx import VcfIdx
from Pipeline import String, Array, File, Tool, ToolOutput, ToolInput


class GatkMutect(Tool):
    tumor = ToolInput("tumor", Bam())
    normal = ToolInput("normal", Bam())
    reference = ToolInput("reference", RefFasta())
    bedFile = ToolInput("bedFile", Bed())
    outputFilename = ToolInput("outputFilename", String())
    dbsnp = ToolInput("dbsnp", Dbsnp())
    cosmic = ToolInput("cosmic", Array(File()))

    vcf = ToolOutput("vcf", VcfIdx())

    @staticmethod
    def tool():
        return "gatk-mutect2"

    @staticmethod
    def base_command():
        return "javac"
