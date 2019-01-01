from abc import ABC

from Pipeline import ToolInput, Array, Filename, ToolOutput, File
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.vcf import VcfIdx
from Pipeline.bioinformatics.tools.gatk3.gatk3toolbase import Gatk3ToolBase


class Gatk3Mutect2Base(Gatk3ToolBase, ABC):

    @staticmethod
    def analysis_type():
        return "MuTect2"

    @staticmethod
    def tool():
        return "gatkmutect2"

    def inputs(self):
        return [
            *super(Gatk3Mutect2Base, self).inputs(),
            *Gatk3Mutect2Base.additional_args,
            ToolInput("inputBam_tumor", BamPair(), position=5, prefix="-I:Tumor",
                      doc="BAM/SAM/CRAM file containing reads, tagged as a 'Tumor'"),
            ToolInput("inputBam_normal", BamPair(), position=6, prefix="-I:Normal",
                      doc="BAM/SAM/CRAM file containing reads, tagged as a 'Normal'"),
            ToolInput("intervals", Bed(), position=7, prefix="-L",
                      doc="One or more genomic intervals over which to operate (previously, .bedFile)"),
            ToolInput("reference", FastaWithDict(), position=8, prefix="-R", doc="Reference sequence file"),
            ToolInput("dbsnp", Array(VcfIdx()), position=10, prefix="--dbsnp"),
            ToolInput("cosmic", Array(VcfIdx()), position=11, prefix="--cosmic",
                      doc="MuTect2 has the ability to use COSMIC data in conjunction with dbSNP to adjust the "
                          "threshold for evidence of a variant in the normal. If a variant is present in dbSNP, "
                          "but not in COSMIC, then more evidence is required from the normal sample to prove the "
                          "variant is not present in germline. This argument supports reference-ordered data (ROD) "
                          "files in the following formats: BCF2, VCF, VCF3"),
            ToolInput("outputFilename", Filename(), position=9, prefix="-o",
                      doc="File to which variants should be written (default: stdout)")
        ]

    def outputs(self):
        return [
            ToolOutput("output", File(), glob="$(inputs.outputFilename)")
        ]

    additional_args = []

    @staticmethod
    def doc():
        return """
    MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping 
    engine of the original MuTect (Cibulskis et al., 2013) with the assembly-based machinery of HaplotypeCaller.
    
    Documentation: https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php
""".strip()


if __name__ == "__main__":
    print(Gatk3Mutect2Base().help())



# class GatkMutect2(CommandTool):
#     inputBam_tumor =
#
#     output = ToolOutput("output", File(), glob='$(inputs.outputfile_name)')
#
#     @staticmethod
#     def tool():
#         return "GatkMutect2"
#
#     @staticmethod
#     def base_command():
#         return ['java']
#
#     @staticmethod
#     def docker():
#         return "broadinstitute/gatk3:3.7-0"
#
#     @staticmethod
#     def doc():
#         return None
#
#     def arguments(self):
#         return [
#             ToolArgument("/usr/GenomeAnalysisTK.jar", position=3, prefix="-jar"),
#             # ToolArgument("/usr/GenomeAnalysisTK.jar", position=3, prefix="-jar"),
#             ToolArgument("MuTect2", position=4, prefix="-T")
#         ]
#
#     java_arg = ToolInput("java_arg", String(optional=True), position=1, default="-Xmx4g")
