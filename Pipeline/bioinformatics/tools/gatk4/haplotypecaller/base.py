from abc import ABC

from Pipeline import String, Int, File, ToolOutput, ToolInput, \
    ToolArgument, Boolean, Double, Array, Filename
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.vcf import VcfIdx, Vcf
from Pipeline.bioinformatics.tools.gatk4.gatk4toolbase import Gatk4ToolBase


class Gatk4HaplotypeCallerBase(Gatk4ToolBase, ABC):
    @classmethod
    def gatk_command(cls):
        return "HaplotypeCaller"

    @staticmethod
    def tool():
        return "GatkHaplotypeCaller"

    def inputs(self):
        return [
            *super(Gatk4HaplotypeCallerBase, self).inputs(),
            *Gatk4HaplotypeCallerBase.optional_args,
            ToolInput("inputRead", BamPair(), doc="BAM/SAM/CRAM file containing reads", prefix="--input"),
            ToolInput("reference", FastaWithDict(), position=5, prefix="--reference", doc="Reference sequence file"),
            ToolInput("outputFilename", Filename(extension=".vcf"), position=8, prefix="--output",
                      doc="File to which variants should be written"),
            ToolInput("dbsnp", VcfIdx(), position=7, prefix="--dbsnp", doc="(Also: -D) A dbSNP VCF file."),
        ]

    def outputs(self):
        return [
            ToolOutput("output", Vcf(), glob='$(inputs.outputFilename)',
                       doc="A raw, unfiltered, highly sensitive callset in VCF format. "
                           "File to which variants should be written"),
        ]

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
    
    Documentation: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php#
""".strip()

    optional_args = [
        ToolInput("activityProfileOut", String(optional=True), prefix="--activity-profile-out",
                  doc="Output the raw activity profile results in IGV format (default: null)"),
        ToolInput("alleles", File(optional=True), prefix="--alleles",
                  doc="(default: null) The set of alleles at which to genotype when --genotyping_mode "
                      "is GENOTYPE_GIVEN_ALLELES"),
        ToolInput("annotateWithNumDiscoveredAlleles", Boolean(optional=True), prefix="--annotate-with-num-discovered-alleles",
                  doc="If provided, we will annotate records with the number of alternate alleles that were "
                      "discovered (but not necessarily genotyped) at a given site"),
        ToolInput("annotation", Array(String(), optional=True), prefix="--annotation",
                  doc="-A: One or more specific annotations to add to variant calls"),
        ToolInput("annotationGroup", Array(String(), optional=True), prefix="--annotation-group",
                  doc="-G	One or more groups of annotations to apply to variant calls"),
        ToolInput("annotationsToExclude", Array(String(), optional=True), prefix="--annotations-to-exclude",
                  doc="-AX	One or more specific annotations to exclude from variant calls"),
        ToolInput("arguments_file", Array(File(), optional=True), prefix="--arguments_file",
                  doc="read one or more arguments files and add them to the command line"),
        ToolInput("assemblyRegionOut", String(optional=True), prefix="--assembly-region-out",
                  doc="(default: null) Output the assembly region to this IGV formatted file. Which annotations to "
                      "exclude from output in the variant calls. Note that this argument has higher priority than "
                      "the -A or -G arguments, so these annotations will be excluded even if they are explicitly "
                      "included with the other options."),
        ToolInput("baseQualityScoreThreshold", Int(optional=True), prefix="--base-quality-score-threshold",
                  doc="(default: 18) Base qualities below this threshold will be reduced to the minimum (6)"),
        ToolInput("cloudIndexPrefetchBuffer", Int(optional=True), prefix="--cloud-index-prefetch-buffer",
                  doc="-CIPB (default: -1) Size of the cloud-only prefetch buffer (in MB; 0 to disable). "
                      "Defaults to cloudPrefetchBuffer if unset."),
        ToolInput("cloudPrefetchBuffer", Int(optional=True), prefix="--cloud-prefetch-buffer",
                  doc="-CPB (default: 40) Size of the cloud-only prefetch buffer (in MB; 0 to disable)."),
        ToolInput("contaminationFractionToFilter", Double(optional=True), prefix="--contamination-fraction-to-filter",
                  doc="-contamination (default: 0.0) Fraction of contamination in sequencing data "
                      "(for all samples) to aggressively remove"),
        ToolInput("correctOverlappingQuality", Boolean(optional=True), prefix="--correct-overlapping-quality",
                  doc="Undocumented option"),
        ToolInput("dbsnp", VcfIdx(optional=True), prefix="--dbsnp", doc="-D (default: null) dbSNP file"),
        ToolInput("disableBamIndexCaching", Boolean(optional=True), prefix="--disable-bam-index-caching",
                  doc="-DBIC. If true, don't cache bam indexes, this will reduce memory requirements but may harm "
                      "performance if many intervals are specified. Caching is automatically disabled if "
                      "there are no intervals specified."),
        # ToolInput("disableSequenceDictionaryValidation", Boolean(optional=True), prefix="--disable-sequence-dictionary-validation",
        #           doc="If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!"),
        ToolInput("founderId", Array(String(), optional=True), prefix="--founder-id",
                  doc="Samples representing the population \"founders\""),
        # ToolInput("gcsMaxRetries", Int(optional=True), prefix="--gcs-max-retries",
        #           doc="-gcs-retries (default: 20) If the GCS bucket channel errors out, "
        #               "how many times it will attempt to re-initiate the connection"),
        # ToolInput("gcsProjectForRequesterPays", String(), prefix="--gcs-project-for-requester-pays",
        #           doc="Project to bill when accessing \"requester pays\" buckets. If unset, these buckets cannot be accessed."),
        ToolInput("genotypingMode", String(optional=True), prefix="--genotyping-mode",
                  doc="(default: DISCOVERY) Specifies how to determine the alternate alleles to use for genotyping. "
                      "The --genotyping-mode argument is an enumerated type (GenotypingOutputMode), which can have one "
                      "of the following values: DISCOVERY (The genotyper will choose the most likely alternate allele) "
                      "or GENOTYPE_GIVEN_ALLELES (Only the alleles passed by the user should be considered)."),
        # ToolInput("graphOutput", DataType(optional=True), prefix="--graph-output", doc="-graph	null	Write debug assembly graph information to this file"),
        ToolInput("heterozygosity", Double(optional=True), prefix="--heterozygosity",
                  doc="(default: 0.001) Heterozygosity value used to compute prior likelihoods for any locus. The "
                      "expected heterozygosity value used to compute prior probability that a locus is non-reference. "
                      "The default priors are for provided for humans: het = 1e-3 which means that the probability "
                      "of N samples being hom-ref at a site is: 1 - sum_i_2N (het / i) Note that heterozygosity as "
                      "used here is the population genetics concept: "
                      "http://en.wikipedia.org/wiki/Zygosity#Heterozygosity_in_population_genetics . "
                      "That is, a hets value of 0.01 implies that two randomly chosen chromosomes from the population "
                      "of organisms would differ from each other (one being A and the other B) at a rate of 1 in 100 bp. "
                      "Note that this quantity has nothing to do with the likelihood of any given sample having a "
                      "heterozygous genotype, which in the GATK is purely determined by the probability of the observed "
                      "data P(D | AB) under the model that there may be a AB het genotype. The posterior probability "
                      "of this AB genotype would use the het prior, but the GATK only uses this posterior probability "
                      "in determining the prob. that a site is polymorphic. So changing the het parameters only "
                      "increases the chance that a site will be called non-reference across all samples, but doesn't "
                      "actually change the output genotype likelihoods at all, as these aren't posterior probabilities "
                      "at all. The quantity that changes whether the GATK considers the possibility of a het genotype "
                      "at all is the ploidy, which determines how many chromosomes each individual in the species carries."),
        ToolInput("heterozygosityStdev", Double(optional=True), prefix="--heterozygosity-stdev",
                  doc="(default 0.01) Standard deviation of heterozygosity for SNP and indel calling."),
        ToolInput("indelHeterozygosity", Double(optional=True), prefix="--indel-heterozygosity",
                  doc="(default: 1.25E-4) Heterozygosity for indel calling. This argument informs the prior "
                      "probability of having an indel at a site. (See heterozygosity)"),
        ToolInput("intervalMergingRule", String(optional=True), prefix="--interval-merging-rule",
                  doc="-imr (default: ALL) Interval merging rule for abutting intervals. By default, the program "
                      "merges abutting intervals (i.e. intervals that are directly side-by-side but do not actually "
                      "overlap) into a single continuous interval. However you can change this behavior if you want "
                      "them to be treated as separate intervals instead. The --interval-merging-rule argument is an "
                      "enumerated type (IntervalMergingRule), which can have one of the following values:"
                      "[ALL, OVERLAPPING]"),
        ToolInput("intervals", Bed(optional=True), prefix="--intervals",
                  doc="-L	One or more genomic intervals over which to operate"),
        ToolInput("maxReadsPerAlignmentStart", Int(optional=True), prefix="--max-reads-per-alignment-start",
                  doc="(default: 50) Maximum number of reads to retain per alignment start position. "
                      "Reads above this threshold will be downsampled. Set to 0 to disable."),
        ToolInput("minBaseQualityScore", Int(optional=True), prefix="--min-base-quality-score",
                  doc="-mbq (default: 10) Minimum base quality required to consider a base for calling"),
        ToolInput("nativePairHmmThreads", Int(optional=True), prefix="--native-pair-hmm-threads",
                  doc="(default: 4) How many threads should a native pairHMM implementation use"),
        ToolInput("nativePairHmmUseDoublePrecision", Boolean(optional=True), prefix="--native-pair-hmm-use-double-precision",
                  doc="use double precision in the native pairHmm. "
                      "This is slower but matches the java implementation better"),
        ToolInput("numReferenceSamplesIfNoCall", Int(optional=True), prefix="--num-reference-samples-if-no-call",
                  doc="(default: 0) Number of hom-ref genotypes to infer at sites not present in a panel. When a "
                      "variant is not seen in any panel, this argument controls whether to infer (and with what "
                      "effective strength) that only reference alleles were observed at that site. "
                      "E.g. \"If not seen in 1000Genomes, treat it as AC=0, AN=2000\"."),
        ToolInput("outputMode", String(optional=True), prefix="--output-mode",
                  doc="(default: EMIT_VARIANTS_ONLY) Specifies which type of calls we should output. The --output-mode "
                      "argument is an enumerated type (OutputMode), which can have one of the following values: "
                      "[EMIT_VARIANTS_ONLY (produces calls only at variant sites), "
                      "EMIT_ALL_CONFIDENT_SITES (produces calls at variant sites and confident reference sites), "
                      "EMIT_ALL_SITES (produces calls at any callable site regardless of confidence; "
                      "this argument is intended only for point mutations (SNPs) in DISCOVERY mode or "
                      "generally when running in GENOTYPE_GIVEN_ALLELES mode; it will by no means produce "
                      "a comprehensive set of indels in DISCOVERY mode)]"),
        ToolInput("pedigree", File(optional=True), prefix="--pedigree",
                  doc="-ped (default: null) Pedigree file for determining the population \"founders\""),
        ToolInput("populationCallset", File(optional=True), prefix="--population-callset",
                  doc="-population (default: null) Callset to use in calculating genotype priors"),
        ToolInput("sampleName", String(optional=True), prefix="--sample-name",
                  doc="-ALIAS (default: null) Name of single sample to use from a multi-sample bam. You can use this "
                      "argument to specify that HC should process a single sample out of a multisample BAM file. "
                      "This is especially useful if your samples are all in the same file but you need to run them "
                      "individually through HC in -ERC GVC mode (which is the recommended usage). "
                      "Note that the name is case-sensitive."),
        ToolInput("samplePloidy", Int(optional=True), prefix="--sample-ploidy",
                  doc="-ploidy (default: 2) Ploidy (number of chromosomes) per sample. "
                      "For pooled data, set to (Number of samples in each pool * Sample Ploidy). "
                      "Sample ploidy - equivalent to number of chromosomes per pool. In pooled "
                      "experiments this should be = # of samples in pool * individual sample ploidy"),
        ToolInput("sitesOnlyVcfOutput", Boolean(optional=True), prefix="--sites-only-vcf-output",
                  doc="(default: false) If true, don't emit genotype fields when writing vcf file output."),
        ToolInput("standardMinConfidenceThresholdForCalling", Double(optional=True),
                  prefix="--standard-min-confidence-threshold-for-calling",
                  doc="-stand-call-conf (default: 10.0) The minimum phred-scaled confidence "
                      "threshold at which variants should be called"),
        ToolInput("useNewQualCalculator", Boolean(optional=True), prefix="--use-new-qual-calculator",
                  doc="-new-qual If provided, we will use the new AF model instead of the so-called exact model")
    ]
