from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline import String, Int, File, CommandTool, ToolOutput, ToolInput, \
    ToolArgument, Boolean, Double, Array, Filename
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.vcfidx import VcfIdx


class GatkHaplotypeCaller(CommandTool):
    reference = ToolInput("reference", FastaWithDict(), position=5, prefix="-R")
    dbsnp = ToolInput("dbsnp", VcfIdx(), position=7, prefix="--dbsnp", doc="latest_dbsnp.vcf set of known indels")
    inputBam_HaplotypeCaller = ToolInput("inputBam_HaplotypeCaller", BamPair(), position=6, prefix="-I",
                                         doc="bam file produced after printReads")
    bedFile = ToolInput("bedFile", Bed(), position=60, prefix="-L")

    outputFilename = ToolInput("outputFilename", String(optional=True), position=8, prefix="-o",
                                           doc="name of the output file from HaplotypeCaller")
    bamOutput = ToolInput("bamOutput", Filename(), position=48, prefix="--bamOutput",
                          doc="File to which assembled haplotypes should be written")

    output = ToolOutput("output", File(), glob='$(inputs.outputFilename)')
    bamOut = ToolOutput("bamOut", File(), glob='$(inputs.bamOutput)')

    @staticmethod
    def tool():
        return "GatkHaplotypeCaller"

    @staticmethod
    def base_command():
        return ['java']

    @staticmethod
    def docker():
        return "broadinstitute/gatk3:3.7-0"

    @staticmethod
    def doc():
        return "None"

    def arguments(self):
        return [ToolArgument("./test/test-files", position=2, prefix="-Djava.io.tmpdir=",
                             separate_value_from_prefix=False),
                ToolArgument("/usr/GenomeAnalysisTK.jar", position=3, prefix="-jar"),
                ToolArgument("HaplotypeCaller", position=4, prefix="-T")]

    max_alternate_alleles = ToolInput("max_alternate_alleles", Int(optional=True), position=25,
                                      prefix="--max_alternate_alleles",
                                      doc="Maximum number of alternate alleles to genotype")
    activeProbabilityThreshold = ToolInput("activeProbabilityThreshold", Double(optional=True), position=58,
                                           prefix="--activeProbabilityThreshold",
                                           doc="Threshold for the probability of a profile state being active.")
    alleles = ToolInput("alleles", Array(String(), optional=True), position=53,
                        doc="The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES")

    stand_emit_conf = ToolInput("stand_emit_conf", Double(optional=True), position=11,
                                prefix="--standard_min_confidence_threshold_for_emitting",
                                doc="The minimum phred-scaled confidence threshold at which variants should be emitted "
                                    "(and filtered with LowQual if less than the calling threshold)")
    kmerSize = ToolInput("kmerSize", Array(Int(), optional=True), position=26,
                         doc="Kmer size to use in the read threading assembler")
    minDanglingBranchLength = ToolInput("minDanglingBranchLength", Int(optional=True), position=21,
                                        prefix="--minDanglingBranchLength",
                                        doc="Minimum length of a dangling branch to attempt recovery")
    bandPassSigma = ToolInput("bandPassSigma", Double(optional=True), position=46, prefix="--consensus",
                              doc="The sigma of the band pass filter Gaussian kernel; "
                                  "if not provided defaults to Walker annotated default")
    maxReadsInRegionPerSample = ToolInput("maxReadsInRegionPerSample", Int(optional=True), position=23,
                                          prefix="--maxReadsInRegionPerSample", doc="Maximum reads in an active region")
    dontIncreaseKmerSizesForCycles = ToolInput("dontIncreaseKmerSizesForCycles", Boolean(optional=True), position=39,
                                               prefix="--dontIncreaseKmerSizesForCycles",
                                               doc="Disable iterating over kmer sizes when graph cycles are detected")
    globalMAPQ = ToolInput("globalMAPQ", Int(optional=True), position=15,
                           prefix="--phredScaledGlobalReadMismappingRate",
                           doc="The global assumed mismapping rate for reads")
    java_arg = ToolInput("java_arg", String(optional=True), position=1, default="-Xmx4g")
    min_base_quality_score = ToolInput("min_base_quality_score", Int(optional=True), position=22,
                                       prefix="--min_base_quality_score",
                                       doc="Minimum base quality required to consider a base for calling")
    excludeAnnotation = ToolInput("excludeAnnotation", Array(String(), optional=True), position=37,
                                  doc="One or more specific annotations to exclude")
    allowNonUniqueKmersInRef = ToolInput("allowNonUniqueKmersInRef", Boolean(optional=True), position=52,
                                         prefix="--allowNonUniqueKmersInRef",
                                         doc="Allow graphs that have non-unique kmers in the reference")
    group = ToolInput("group", Array(String(), optional=True), position=32, doc="Input prior for calls")
    pcr_indel_model = ToolInput("pcr_indel_model", String(optional=True), position=16, prefix="--pcr_indel_model",
                                doc="The PCR indel model to use")
    stand_call_conf = ToolInput("stand_call_conf", Double(optional=True), position=12,
                                prefix="--standard_min_confidence_threshold_for_calling",
                                doc="The minimum phred-scaled confidence threshold at which variants should be called")
    activeRegionExtension = ToolInput("activeRegionExtension", Int(optional=True), position=57,
                                      prefix="--activeRegionExtension",
                                      doc="The active region extension; "
                                          "if not provided defaults to Walker annotated default")
    activeRegionOut = ToolInput("activeRegionOut", File(optional=True), position=55, prefix="--activeRegionOut",
                                doc="Output the active region to this IGV formatted file")
    useAllelesTrigger = ToolInput("useAllelesTrigger", Boolean(optional=True), position=10,
                                  prefix="--useAllelesTrigger",
                                  doc="Use additional trigger on variants found in an external alleles file")
    forceActive = ToolInput("forceActive", Boolean(optional=True), position=36, prefix="--forceActive",
                            doc="If provided, all bases will be tagged as active")
    sample_name = ToolInput("sample_name", String(optional=True), position=14, prefix="--sample_name",
                            doc="Use additional trigger on variants found in an external alleles file")
    useFilteredReadsForAnnotations = ToolInput("useFilteredReadsForAnnotations", Boolean(optional=True), position=9,
                                               prefix="--useFilteredReadsForAnnotations",
                                               doc="Use the contamination-filtered read maps for the"
                                                   " purposes of annotating variants")
    disableOptimizations = ToolInput("disableOptimizations", Boolean(optional=True), position=41,
                                     prefix="--disableOptimizations",
                                     doc="Dont skip calculations in ActiveRegions with no variants")
    minPruning = ToolInput("minPruning", Int(optional=True), position=20, prefix="--minPruning",
                           doc="Minimum support to not prune paths in the graph")
    activeRegionMaxSize = ToolInput("activeRegionMaxSize", Int(optional=True), position=56,
                                    prefix="--activeRegionMaxSize",
                                    doc="The active region maximum size; "
                                        "if not provided defaults to Walker annotated default")
    output_mode = ToolInput("output_mode", String(optional=True), position=17, prefix="--output_mode",
                            doc="The PCR indel model to use")
    annotateNDA = ToolInput("annotateNDA", Boolean(optional=True), position=50, prefix="--annotateNDA",
                            doc="If provided, we will annotate records with the number of alternate alleles that "
                                "were discovered (but not necessarily genotyped) at a given site")
    ERCIS = ToolInput("ERCIS", Int(optional=True), position=28, prefix="--indelSizeToEliminateInRefModel",
                      doc="The size of an indel to check for in the reference model")
    GVCFGQBands = ToolInput("GVCFGQBands", Array(Int(), optional=True), position=31, doc="Input prior for calls")
    allSitePLs = ToolInput("allSitePLs", Boolean(optional=True), position=51, prefix="--allSitePLs",
                           doc="Annotate all sites with PLs")
    numPruningSamples = ToolInput("numPruningSamples", Int(optional=True), position=18, prefix="--numPruningSamples",
                                  doc="Number of samples that must pass the minPruning threshold")
    gcpHMM = ToolInput("gcpHMM", Int(optional=True), position=35, prefix="--gcpHMM",
                       doc="Flat gap continuation penalty for use in the Pair HMM")
    contamination = ToolInput("contamination", File(optional=True), position=43,
                              prefix="--contamination_fraction_to_filter",
                              doc="Tab-separated File containing fraction of contamination in sequencing data "
                                  "(per sample) to aggressively remove. "
                                  "Format should be (Contamination is double) per line; No header.")
    graphOutput = ToolInput("graphOutput", File(optional=True), position=33, prefix="--graphOutput",
                            doc="Write debug assembly graph information to this file")
    dontTrimActiveRegions = ToolInput("dontTrimActiveRegions", Boolean(optional=True), position=39,
                                      prefix="--dontTrimActiveRegions",
                                      doc="If specified, we will not trim down the active region from the full "
                                          "region (active + extension) to just the active interval for genotyping")
    annotation = ToolInput("annotation", Array(String(), optional=True), position=49,
                           doc="One or more specific annotations to apply to variant calls")
    bamWriterType = ToolInput("bamWriterType", String(optional=True), position=47, prefix="--bamWriterType",
                              doc="Which haplotypes should be written to the BAM.")
    genotyping_mode = ToolInput("genotyping_mode", String(optional=True), position=34, prefix="--genotyping_mode",
                                doc="The --genotyping_mode argument is an enumerated type (GenotypingOutputMode), "
                                    "which can have one of the following values")
    activityProfileOut = ToolInput("activityProfileOut", File(optional=True), position=54,
                                   prefix="--activityProfileOut",
                                   doc="Output the raw activity profile results in IGV format")
    input_prior = ToolInput("input_prior", Array(Double(), optional=True), position=27, doc="Input prior for calls")
    indel_heterozygosity = ToolInput("indel_heterozygosity", Double(optional=True), position=29,
                                     prefix="--indel_heterozygosity", doc="Heterozygosity for indel calling")
    emitRefConfidenceDBSN = ToolInput("emitRefConfidenceDBSN", String(optional=True), position=38,
                                      prefix="--emitRefConfidenceDBSN",
                                      doc="Mode for emitting reference confidence scores")
    consensus = ToolInput("consensus", Boolean(optional=True), position=44, prefix="--consensus",
                          doc="Print out very verbose debug information about each triggering active region")
    heterozygosity = ToolInput("heterozygosity", Double(optional=True), position=30, prefix="--heterozygosity",
                               doc="Heterozygosity for indel calling")
    minReadsPerAlignmentStart = ToolInput("minReadsPerAlignmentStart", Int(optional=True), position=19,
                                          prefix="--minReadsPerAlignmentStart",
                                          doc="Minimum number of reads sharing the same alignment start for each "
                                              "genomic location in an active region")
    sample_ploidy = ToolInput("sample_ploidy", Int(optional=True), position=13, prefix="--sample_ploidy",
                              doc="Use additional trigger on variants found in an external alleles file")
    debug = ToolInput("debug", Boolean(optional=True), position=42, prefix="--debug",
                      doc="Print out very verbose debug information about each triggering active region")
    doNotRunPhysicalPhasing = ToolInput("doNotRunPhysicalPhasing", Boolean(optional=True), position=40,
                                        prefix="--doNotRunPhysicalPhasing",
                                        doc="As of GATK 3.3, HaplotypeCaller outputs physical (read-based) information "
                                            "(see version 3.3 release notes and documentation for details). "
                                            "This argument disables that behavior.")
    comp = ToolInput("comp", Array(String(), optional=True), position=45,
                     doc="comp binds reference ordered data. "
                         "This argument supports ROD files of the following types BCF2, VCF, VCF3")
    maxNumHaplotypesInPopulation = ToolInput("maxNumHaplotypesInPopulation", Int(optional=True), position=24,
                                             prefix="--maxNumHaplotypesInPopulation",
                                             doc="Maximum number of haplotypes to consider for your population")
    threads = ToolInput("threads", Int(optional=True), position=56, prefix="-nct")
    emitRefConfidence = ToolInput("emitRefConfidence", String(optional=True), position=61, prefix="--emitRefConfidence",
                                  default="NONE")


if __name__ == "__main__":
    import yaml

    print(GatkHaplotypeCaller().help())
    # print(yaml.dump(GatkHaplotypeCaller().cwl(), default_flow_style=False))
