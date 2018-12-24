from abc import ABC

from Pipeline import ToolInput, Array, Filename, ToolArgument, ToolOutput
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.vcf import VcfIdx
from Pipeline.bioinformatics.tools.gatk.gatk4base import Gatk4Base


class GatkMutect2Base(Gatk4Base, ABC):
    @staticmethod
    def tool():
        return "gatkmutect2"

    def inputs(self):
        return [
            *super(GatkMutect2Base, self).inputs(),
            *GatkMutect2Base.additional_args,
            ToolInput("inputBam_tumor", BamPair(), position=5, prefix="-I:Tumor",
                      doc="BAM/SAM/CRAM file containing reads, tagged as a 'Tumor'"),
            ToolInput("inputBam_normal", BamPair(), position=6, prefix="-I:Normal",
                      doc="BAM/SAM/CRAM file containing reads, tagged as a 'Normal'"),
            ToolInput("intervals", Bed(), position=7, prefix="-L",
                      doc="One or more genomic intervals over which to operate (previously, .bedFile)"),
            ToolInput("reference", FastaWithDict(), position=8, prefix="-R", doc="Reference sequence file"),
            ToolInput("dbsnp", Array(VcfIdx()), position=10, prefix="--dbsnp"),
            ToolInput("cosmic", Array(VcfIdx()), position=11, prefix="--cosmic"),   # idk about this?
            ToolInput("outputFilename", Filename(), position=9, prefix="-o")
        ]

    def outputs(self):
        return [
            ToolOutput("output", File(), glob="")
        ]

    def arguments(self):
        return [
            *super(GatkMutect2Base, self).arguments(),
            ToolArgument("MuTect2", position=4, prefix="-T")
        ]

    additional_args = []

    @staticmethod
    def doc():
        return """
    Call somatic short variants via local assembly of haplotypes. Short variants include single nucleotide (SNV) 
    and insertion and deletion (indel) variants. The caller combines the DREAM challenge-winning somatic 
    genotyping engine of the original MuTect (Cibulskis et al., 2013) with the assembly-based machinery of HaplotypeCaller.

    This tool is featured in the Somatic Short Mutation calling Best Practice Workflow. See Tutorial#11136 
    for a step-by-step description of the workflow and Article#11127 for an overview of what traditional 
    somatic calling entails. For the latest pipeline scripts, see the Mutect2 WDL scripts directory. 
    Although we present the tool for somatic calling, it may apply to other contexts, 
    such as mitochondrial variant calling.
    
    Documentation: https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.10.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php
""".strip()

if __name__ == "__main__":
    print(GatkMutect2Base().help())



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
