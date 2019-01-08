from abc import ABC

from Pipeline import ToolArgument, ToolInput, Filename, ToolOutput, File, Array, String, Boolean, Int, Float, Directory
from Pipeline.bioinformatics.data_types.vcf import VcfIdx, TabixIdx
from Pipeline.bioinformatics.tools.gatk4.gatk4toolbase import Gatk4ToolBase


class Gatk4GenotypeConcordanceBase(Gatk4ToolBase, ABC):
    @classmethod
    def gatk_command(cls):
        return "GenotypeConcordance"

    @staticmethod
    def tool():
        return "Gatk4GenotypeConcordance"

    def inputs(self):
        return [
            ToolInput("callVCF", TabixIdx(), prefix="--CALL_VCF", doc="The VCF containing the call sample"),
            ToolInput("truthVCF", VcfIdx(), prefix="--TRUTH_VCF", doc="The VCF containing the truth sample"),
            ToolInput("outputBasename", Filename(), prefix="--OUTPUT",
                      doc="Basename for the three metrics files that are to be written. Resulting files will be:"
                          "(1) .genotype_concordance_summary_metrics, "
                          "(2) .genotype_concordance_detail_metrics, "
                          "(3) .genotype_concordance_contingency_metrics."),
            *super(Gatk4GenotypeConcordanceBase, self).inputs(),
            *self.additional_args
        ]

    def outputs(self):
        return [
            ToolOutput("summaryMetrics", File(), glob="*.genotype_concordance_summary_metrics"),
            ToolOutput("detailMetrics", File(), glob="*.genotype_concordance_detail_metrics"),
            ToolOutput("contingencyMetrics", File(), glob="*.genotype_concordance_contingency_metrics"),
            # ToolOutput("vcf", VcfIdx(optional=True), glob="*.vcf")
        ]

    @staticmethod
    def doc():

        return """
    GenotypeConcordance (Picard)
    Calculates the concordance between genotype data of one samples in each of two VCFs - one being 
    considered the truth (or reference) the other being the call. The concordance is broken into 
    separate results sections for SNPs and indels. Statistics are reported in three different files.
    
    Summary
    Calculates the concordance between genotype data of one samples in each of two VCFs - one being 
    considered the truth (or reference) the other being the call. The concordance is broken into 
    separate results sections for SNPs and indels. Summary and detailed statistics are reported.
    
    Details
    This tool evaluates the concordance between genotype calls for a sample in different callsets
    where one is being considered as the "truth" (aka standard, or reference) and the other as the 
    "call" that is being evaluated for accuracy. The Comparison can be restricted to a confidence 
    interval which is typically used in order to enable proper assessment of False Positives and 
    the False-Positive Rate (FPR).
     
    Output Metrics:
    Output metrics consists of GenotypeConcordanceContingencyMetrics, GenotypeConcordanceSummaryMetrics, 
    and GenotypeConcordanceDetailMetrics. For each set of metrics, the data is broken into separate 
    sections for SNPs and INDELs. Note that only SNP and INDEL variants are considered, MNP, Symbolic, 
    and Mixed classes of variants are not included.
    
    GenotypeConcordanceContingencyMetrics enumerate the constituents of each contingent in a callset 
    including true-positive (TP), true-negative (TN), false-positive (FP), and false-negative (FN) calls.
    GenotypeConcordanceDetailMetrics include the numbers of SNPs and INDELs for each contingent genotype 
    as well as the number of validated genotypes.
    
    GenotypeConcordanceSummaryMetrics provide specific details for the variant caller performance 
    on a callset including values for sensitivity, specificity, and positive predictive values.
    
    
    Useful definitions applicable to alleles and genotypes:
        - Truthset - A callset (typically in VCF format) containing variant calls and genotypes that have been 
            cross-validated with multiple technologies e.g. Genome In A Bottle Consortium (GIAB) (https://sites.stanford.edu/abms/giab)
        - TP - True-positives are variant sites that match against the truth-set
        - FP - False-positives are reference sites miscalled as variant
        - FN - False-negatives are variant sites miscalled as reference
        - TN - True-negatives are correctly called as reference
        - Validated genotypes - are TP sites where the exact genotype (HET or HOM-VAR) appears in the truth-set
    
    VCF Output:
        - The concordance state will be stored in the CONC_ST tag in the INFO field
        - The truth sample name will be \"truth\" and call sample name will be \"call\"  
        
    Documentation: https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_vcf_GenotypeConcordance.php
    """.strip()

    additional_args = [
        ToolInput("argumentsFile", Array(File(), optional=True), prefix="--arguments_file", position=10,
                  doc="read one or more arguments files and add them to the command line"),
        ToolInput("callSample", String(optional=True), prefix="--CALL_SAMPLE", position=10,
                  doc="The name of the call sample within the call VCF. Not required if only one sample exists."),

        ToolInput("ignoreFilterStatus", Boolean(optional=True), prefix="--IGNORE_FILTER_STATUS",
                  doc="Default is false. If true, filter status of sites will be ignored so that we "
                      "include filtered sites when calculating genotype concordance."),
        ToolInput("intersectIntervals", Boolean(optional=True), prefix="--INTERSECT_INTERVALS",
                  doc="If true, multiple interval lists will be intersected. If false multiple lists will be unioned."),
        ToolInput("intervals", Array(VcfIdx(), optional=True), prefix="--INTERVALS",
                  doc="One or more interval list files that will be used to limit the genotype concordance. "
                      "Note - if intervals are specified, the VCF files must be indexed."),
        ToolInput("minDP", Float(optional=True), prefix="--MIN_DP",
                  doc="Genotypes below this depth will have genotypes classified as LowDp."),
        ToolInput("minGQ", Float(optional=True), prefix="--MIN_GQ",
                  doc="Genotypes below this genotype quality will have genotypes classified as LowGq."),
        ToolInput("treatMissingSitesAsHomeRef", Boolean(optional=True), prefix="--MISSING_SITES_HOM_REF",
                  doc="Default is false, which follows the GA4GH Scheme. If true, missing sites in the truth "
                      "set will be treated as HOM_REF sites and sites missing in both the truth and call sets "
                      "will be true negatives. Useful when hom ref sites are left out of the truth set. "
                      "This flag can only be used with a high confidence interval list."),
        ToolInput("outputAllRows", Boolean(optional=True), prefix="--OUTPUT_ALL_ROWS",
                  doc="If true, output all rows in detailed statistics even when count == 0. When false only "
                      "output rows with non-zero counts."),
        ToolInput("outputVcf", Boolean(optional=True), prefix="--OUTPUT_VCF",
                  doc="Output a VCF annotated with concordance information."),
        ToolInput("truthSample", String(optional=True), prefix="--TRUTH_SAMPLE",
                  doc="The name of the truth sample within the truth VCF. Not required if only one sample exists."),
        ToolInput("useVcfIndex", Boolean(optional=True), prefix="--USE_VCF_INDEX",
                  doc="If true, use the VCF index, else iterate over the entire VCF"),

        ToolInput("compressionLevel", Int(optional=True), prefix="--COMPRESSION_LEVEL", position=11,
                  doc="Compression level for all compressed files created (e.g. BAM and GELI)."),
        ToolInput("createIndex", Boolean(optional=True), prefix="--CREATE_INDEX", position=11,
                  doc="Whether to create a BAM index when writing a coordinate-sorted BAM file."),
        ToolInput("createMd5File", Boolean(optional=True), prefix="--CREATE_MD5_FILE", position=11,
                  doc="Whether to create an MD5 digest for any BAM or FASTQ files created."),
        ToolInput("maxRecordsInRam", Int(optional=True), prefix="--MAX_RECORDS_IN_RAM", position=11,
                  doc="When writing SAM files that need to be sorted, this will specify the number of "
                      "records stored in RAM before spilling to disk. Increasing this number reduces "
                      "the number of file handles needed to sort a SAM file, and increases the amount of RAM needed."),
        ToolInput("quiet", Boolean(optional=True), prefix="--QUIET", position=11,
                  doc="Whether to suppress job-summary info on System.err."),
        ToolInput("reference", File(optional=True), prefix="--REFERENCE=SEQUENCE", position=11,
                  doc="Reference sequence file."),
        ToolInput("tmpDir", Directory(optional=True), prefix="--TMP_DIR", position=11,
                  doc="Undocumented option"),
        ToolInput("useJdkDeflater", Boolean(optional=True), prefix="--use_jdk_deflater", position=11,
                  doc="Whether to use the JdkDeflater (as opposed to IntelDeflater)"),
        ToolInput("useJdkInflater", Boolean(optional=True), prefix="--use_jdk_inflater", position=11,
                  doc="Whether to use the JdkInflater (as opposed to IntelInflater)"),
        ToolInput("validationStringency", String(optional=True), prefix="--VALIDATION_STRINGENCY", position=11,
                  doc="Validation stringency for all SAM files read by this program. Setting stringency to SILENT "
                      "can improve performance when processing a BAM file in which variable-length data "
                      "(read, qualities, tags) do not otherwise need to be decoded."
                      "The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), "
                      "which can have one of the following values: [STRICT, LENIENT, SILENT]"),
        ToolInput("verbosity", String(optional=True), prefix="--verbosity", position=11,
                  doc="The --verbosity argument is an enumerated type (LogLevel), which can have "
                      "one of the following values: [ERROR, WARNING, INFO, DEBUG]")
    ]
