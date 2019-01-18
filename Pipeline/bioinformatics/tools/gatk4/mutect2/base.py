from abc import ABC

from Pipeline import ToolInput, Array, Filename, ToolArgument, ToolOutput, File, String, Float
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.fasta import FastaWithDict
from Pipeline.bioinformatics.data_types.vcf import VcfIdx, Vcf
from Pipeline.bioinformatics.tools.gatk4.gatk4toolbase import Gatk4ToolBase


class Gatk4Mutect2Base(Gatk4ToolBase, ABC):
    @classmethod
    def gatk_command(cls):
        return "Mutect2"

    @staticmethod
    def tool():
        return "gatkmutect2"

    @staticmethod
    def requirements():
        import cwlgen.cwlgen as cwl
        return [
            cwl.ResourceRequirement(ram_min="64000")
        ]

    @staticmethod
    def tumor_normal_inputs():
        return [
            ToolInput("tumor", BamPair(), position=5, prefix="-I", doc="BAM/SAM/CRAM file containing reads"),
            ToolInput("tumorName", String(), position=6, prefix="-tumor",
                      doc="BAM sample name of tumor. May be URL-encoded as output by GetSampleName with -encode."),
            ToolInput("normal", BamPair(), position=5, prefix="-I", doc="BAM/SAM/CRAM file containing reads"),
            ToolInput("normalName", String(), position=6, prefix="-normal",
                      doc="BAM sample name of normal. May be URL-encoded as output by GetSampleName with -encode."),
        ]

    def inputs(self):
        return [
            *super(Gatk4Mutect2Base, self).inputs(),
            *Gatk4Mutect2Base.additional_args,
            *Gatk4Mutect2Base.tumor_normal_inputs(),
            ToolInput("intervals", Bed(optional=True), position=7, prefix="-L",
                      doc="One or more genomic intervals over which to operate"),
            ToolInput("reference", FastaWithDict(), position=8, prefix="-R", doc="Reference sequence file"),
            ToolInput("outputFilename", Filename(extension=".vcf.gz"), position=20, prefix="-O"),
            ToolInput("germlineResource", VcfIdx(optional=True), position=10, prefix="--germline-resource"),
            ToolInput("afOfAllelesNotInResource", Float(optional=True), position=11,
                      prefix="--af-of-alleles-not-in-resource",
                      doc="Population allele fraction assigned to alleles not found in germline resource. "
                          "Please see docs/mutect/mutect2.pdf fora derivation of the default value."),
            ToolInput("panelOfNormals", VcfIdx(optional=True), position=10, prefix="--panel-of-normals",
                      doc="A panel of normals can be a useful (optional) input to help filter out "
                          "commonly seen sequencing noise that may appear as low allele-fraction somatic variants.")
        ]

    def outputs(self):
        return [
            # Todo: Determine type of Gatk4Mutect2 output (.vcf.gz?)
            ToolOutput("output", Vcf(), glob="$(inputs.outputFilename)", doc="To determine type")
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

# if __name__ == "__main__":
#     print(Gatk4Mutect2Base().help())



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
