from abc import ABC

from Pipeline import String, Int, File, ToolOutput, ToolInput, \
    ToolArgument, Boolean, Double, Array, Filename
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.vcf import VcfIdx, Vcf
from Pipeline.bioinformatics.tools.gatk4.gatk4toolbase import Gatk4ToolBase


class Gatk4HaplotypeCallerBase(Gatk4ToolBase, ABC):

    @staticmethod
    def tool():
        return "GatkHaplotypeCaller"

    def inputs(self):
        return [
            *super(Gatk4HaplotypeCallerBase, self).inputs(),
            *Gatk4HaplotypeCallerBase.optional_args,
            ToolInput("inputRead", Bam(), doc="BAM/SAM/CRAM file containing reads", prefix="--input"),
            ToolInput("reference", FastaWithDict(), position=5, prefix="-R", doc="Reference sequence file"),
            ToolInput("outputFilename", String(optional=True), position=8, prefix="-O",
                      doc="File to which variants should be written"),
            ToolInput("bamOutput", Filename(), position=48, prefix="--bamout",
                      doc="File to which assembled haplotypes should be written (prefix previously --bam-output)"),
            ToolInput("dbsnp", VcfIdx(), position=7, prefix="--dbsnp", doc="A dbSNP VCF file."),
            ToolInput("intervals", Bed(), position=60, prefix="-L",
                      doc="(Previously: .bedFile) One or more genomic intervals over which to operate")
        ]

    def outputs(self):
        return [
            ToolOutput("output", Vcf(), glob='$(inputs.outputFilename)',
                       doc="""
    Either a VCF or GVCF file with raw, unfiltered SNP and indel calls. Regular VCFs must be filtered 
    either by variant recalibration (Best Practice) or hard-filtering before use in downstream analyses. 
    If using the GVCF workflow, the output is a GVCF file that must first be run through GenotypeGVCFs 
    and then filtering before further analysis""".strip()),

            ToolOutput("bamOut", Bam(), glob='$(inputs.bamOutput)',
                       doc="File to which assembled haplotypes should be written")
        ]


    # def base_command(self):
    #     base = super(GatkHaplotypeCallerBase, self).base_command()
    #     if isinstance(base, list): base = base[0]
    #     return [base, "haplotypecaller"]


    @staticmethod
    def doc():
        return """
    Call germline SNPs and indels via local re-assembly of haplotypes
    
    The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes 
    in an active region. In other words, whenever the program encounters a region showing signs of variation, it 
    discards the existing mapping information and completely reassembles the reads in that region. This allows the 
    HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when 
    they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at 
    calling indels than position-based callers like UnifiedGenotyper.
    
    In the GVCF workflow used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to 
    generate an intermediate GVCF (not to be used in final analysis), which can then be used in GenotypeGVCFs for joint 
    genotyping of multiple samples in a very efficient way. The GVCF workflow enables rapid incremental processing of 
    samples as they roll off the sequencer, as well as scaling to very large cohort sizes (e.g. the 92K exomes of ExAC).
    
    In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. 
    Note however that the algorithms used to calculate variant likelihoods is not well suited to extreme allele 
    frequencies (relative to ploidy) so its use is not recommended for somatic (cancer) variant discovery. 
    For that purpose, use Mutect2 instead.
    
    Finally, HaplotypeCaller is also able to correctly handle the splice junctions that make RNAseq a challenge 
    for most variant callers, on the condition that the input read data has previously been processed according 
    to our recommendations as documented (https://software.broadinstitute.org/gatk/documentation/article?id=4067).
    
    Documentation: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php#--intervals
""".strip()

    def arguments(self):
        return [*super(Gatk4HaplotypeCallerBase, self).arguments(),
                ToolArgument("HaplotypeCaller", position=4, prefix="-T")]

    optional_args = [
        ToolInput("max_alternate_alleles", Int(optional=True), position=25,
                  prefix="--max_alternate_alleles",
                  doc="Maximum number of alternate alleles to genotype"),
        ToolInput("activeProbabilityThreshold", Double(optional=True), position=58,
                  prefix="--activeProbabilityThreshold",
                  doc="Threshold for the probability of a profile state being active."),
        ToolInput("alleles", Array(String(), optional=True), position=53,
                  doc="The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES"),

        ToolInput("stand_emit_conf", Double(optional=True), position=11,
                  prefix="--standard_min_confidence_threshold_for_emitting",
                  doc="The minimum phred-scaled confidence threshold at which variants should be emitted "
                      "(and filtered with LowQual if less than the calling threshold)"),
        ToolInput("kmerSize", Array(Int(), optional=True), position=26,
                  doc="Kmer size to use in the read threading assembler"),
        ToolInput("minDanglingBranchLength", Int(optional=True), position=21,
                  prefix="--minDanglingBranchLength",
                  doc="Minimum length of a dangling branch to attempt recovery"),
        ToolInput("bandPassSigma", Double(optional=True), position=46, prefix="--consensus",
                  doc="The sigma of the band pass filter Gaussian kernel; "
                      "if not provided defaults to Walker annotated default"),
        ToolInput("maxReadsInRegionPerSample", Int(optional=True), position=23,
                  prefix="--maxReadsInRegionPerSample", doc="Maximum reads in an active region"),
        ToolInput("dontIncreaseKmerSizesForCycles", Boolean(optional=True), position=39,
                  prefix="--dontIncreaseKmerSizesForCycles",
                  doc="Disable iterating over kmer sizes when graph cycles are detected"),
        ToolInput("globalMAPQ", Int(optional=True), position=15,
                  prefix="--phredScaledGlobalReadMismappingRate",
                  doc="The global assumed mismapping rate for reads"),
        ToolInput("java_arg", String(optional=True), position=1, default="-Xmx4g"),
        ToolInput("min_base_quality_score", Int(optional=True), position=22,
                  prefix="--min_base_quality_score",
                  doc="Minimum base quality required to consider a base for calling"),
        ToolInput("excludeAnnotation", Array(String(), optional=True), position=37,
                  doc="One or more specific annotations to exclude"),
        ToolInput("allowNonUniqueKmersInRef", Boolean(optional=True), position=52,
                  prefix="--allowNonUniqueKmersInRef",
                  doc="Allow graphs that have non-unique kmers in the reference"),
        ToolInput("group", Array(String(), optional=True), position=32, doc="Input prior for calls"),
        ToolInput("pcr_indel_model", String(optional=True), position=16, prefix="--pcr_indel_model",
                  doc="The PCR indel model to use"),
        ToolInput("stand_call_conf", Double(optional=True), position=12,
                  prefix="--standard_min_confidence_threshold_for_calling",
                  doc="The minimum phred-scaled confidence threshold at which variants should be called"),
        ToolInput("activeRegionExtension", Int(optional=True), position=57,
                  prefix="--activeRegionExtension",
                  doc="The active region extension; "
                      "if not provided defaults to Walker annotated default"),
        ToolInput("activeRegionOut", File(optional=True), position=55, prefix="--activeRegionOut",
                  doc="Output the active region to this IGV formatted file"),
        ToolInput("useAllelesTrigger", Boolean(optional=True), position=10,
                  prefix="--useAllelesTrigger",
                  doc="Use additional trigger on variants found in an external alleles file"),
        ToolInput("forceActive", Boolean(optional=True), position=36, prefix="--forceActive",
                  doc="If provided, all bases will be tagged as active"),
        ToolInput("sample_name", String(optional=True), position=14, prefix="--sample_name",
                  doc="Use additional trigger on variants found in an external alleles file"),
        ToolInput("useFilteredReadsForAnnotations", Boolean(optional=True), position=9,
                  prefix="--useFilteredReadsForAnnotations",
                  doc="Use the contamination-filtered read maps for the"
                      " purposes of annotating variants"),
        ToolInput("disableOptimizations", Boolean(optional=True), position=41,
                  prefix="--disableOptimizations",
                  doc="Dont skip calculations in ActiveRegions with no variants"),
        ToolInput("minPruning", Int(optional=True), position=20, prefix="--minPruning",
                  doc="Minimum support to not prune paths in the graph"),
        ToolInput("activeRegionMaxSize", Int(optional=True), position=56,
                  prefix="--activeRegionMaxSize",
                  doc="The active region maximum size; "
                      "if not provided defaults to Walker annotated default"),
        ToolInput("output_mode", String(optional=True), position=17, prefix="--output_mode",
                  doc="The PCR indel model to use"),
        ToolInput("annotateNDA", Boolean(optional=True), position=50, prefix="--annotateNDA",
                  doc="If provided, we will annotate records with the number of alternate alleles that "
                      "were discovered (but not necessarily genotyped) at a given site"),
        ToolInput("ERCIS", Int(optional=True), position=28, prefix="--indelSizeToEliminateInRefModel",
                  doc="The size of an indel to check for in the reference model"),
        ToolInput("GVCFGQBands", Array(Int(), optional=True), position=31, doc="Input prior for calls"),
        ToolInput("allSitePLs", Boolean(optional=True), position=51, prefix="--allSitePLs",
                  doc="Annotate all sites with PLs"),
        ToolInput("numPruningSamples", Int(optional=True), position=18, prefix="--numPruningSamples",
                  doc="Number of samples that must pass the minPruning threshold"),
        ToolInput("gcpHMM", Int(optional=True), position=35, prefix="--gcpHMM",
                  doc="Flat gap continuation penalty for use in the Pair HMM"),
        ToolInput("contamination", File(optional=True), position=43,
                  prefix="--contamination_fraction_to_filter",
                  doc="Tab-separated File containing fraction of contamination in sequencing data "
                      "(per sample) to aggressively remove. "
                      "Format should be (Contamination is double) per line; No header."),
        ToolInput("graphOutput", File(optional=True), position=33, prefix="--graphOutput",
                  doc="Write debug assembly graph information to this file"),
        ToolInput("dontTrimActiveRegions", Boolean(optional=True), position=39,
                  prefix="--dontTrimActiveRegions",
                  doc="If specified, we will not trim down the active region from the full "
                      "region (active + extension) to just the active interval for genotyping"),
        ToolInput("annotation", Array(String(), optional=True), position=49,
                  doc="One or more specific annotations to apply to variant calls"),
        ToolInput("bamWriterType", String(optional=True), position=47, prefix="--bamWriterType",
                  doc="Which haplotypes should be written to the BAM."),
        ToolInput("genotyping_mode", String(optional=True), position=34, prefix="--genotyping_mode",
                  doc="The --genotyping_mode argument is an enumerated type (GenotypingOutputMode), "
                      "which can have one of the following values"),
        ToolInput("activityProfileOut", File(optional=True), position=54,
                  prefix="--activityProfileOut",
                  doc="Output the raw activity profile results in IGV format"),
        ToolInput("input_prior", Array(Double(), optional=True), position=27, doc="Input prior for calls"),
        ToolInput("indel_heterozygosity", Double(optional=True), position=29,
                  prefix="--indel_heterozygosity", doc="Heterozygosity for indel calling"),
        ToolInput("emitRefConfidenceDBSN", String(optional=True), position=38,
                  prefix="--emitRefConfidenceDBSN",
                  doc="Mode for emitting reference confidence scores"),
        ToolInput("consensus", Boolean(optional=True), position=44, prefix="--consensus",
                  doc="Print out very verbose debug information about each triggering active region"),
        ToolInput("heterozygosity", Double(optional=True), position=30, prefix="--heterozygosity",
                  doc="Heterozygosity for indel calling"),
        ToolInput("minReadsPerAlignmentStart", Int(optional=True), position=19,
                  prefix="--minReadsPerAlignmentStart",
                  doc="Minimum number of reads sharing the same alignment start for each "
                      "genomic location in an active region"),
        ToolInput("sample_ploidy", Int(optional=True), position=13, prefix="--sample_ploidy",
                  doc="Use additional trigger on variants found in an external alleles file"),
        ToolInput("debug", Boolean(optional=True), position=42, prefix="--debug",
                  doc="Print out very verbose debug information about each triggering active region"),
        ToolInput("doNotRunPhysicalPhasing", Boolean(optional=True), position=40,
                  prefix="--doNotRunPhysicalPhasing",
                  doc="As of GATK 3.3, HaplotypeCaller outputs physical (read-based) information "
                      "(see version 3.3 release notes and documentation for details). "
                      "This argument disables that behavior."),
        ToolInput("comp", Array(String(), optional=True), position=45,
                  doc="comp binds reference ordered data. "
                      "This argument supports ROD files of the following types BCF2, VCF, VCF3"),
        ToolInput("maxNumHaplotypesInPopulation", Int(optional=True), position=24,
                  prefix="--maxNumHaplotypesInPopulation",
                  doc="Maximum number of haplotypes to consider for your population"),
        ToolInput("threads", Int(optional=True), position=56, prefix="-nct"),
        ToolInput("emitRefConfidence", String(optional=True), position=61, prefix="--emitRefConfidence",
                  default="NONE")
    ]


if __name__ == "__main__":
    print(Gatk4HaplotypeCallerBase().help())
    # print(yaml.dump(GatkHaplotypeCaller().cwl(), default_flow_style=False))
