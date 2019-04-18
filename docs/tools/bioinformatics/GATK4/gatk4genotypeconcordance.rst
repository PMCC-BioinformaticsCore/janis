
GATK4: Genotype Concordance
======================================================
Tool identifier: ``gatk4genotypeconcordance``

Tool path: ``from janis_bioinformatics.tools.gatk4 import Gatk4GenotypeConcordance_4_0``

Documentation
-------------

Docker
******
``broadinstitute/gatk:4.0.12.0``

URL
******
`https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_vcf_GenotypeConcordance.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_vcf_GenotypeConcordance.php>`_

Docstring
*********
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
    - The truth sample name will be "truth" and call sample name will be "call"

Outputs
-------
==================  ======  ===============
name                type    documentation
==================  ======  ===============
summaryMetrics      File
detailMetrics       File
contingencyMetrics  File
==================  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

========  ==========  ===========  ==========  ===================================
name      type        prefix       position    documentation
========  ==========  ===========  ==========  ===================================
callVCF   vcf-gz-tbi  --CALL_VCF               The VCF containing the call sample
truthVCF  VCFIDX      --TRUTH_VCF              The VCF containing the truth sample
========  ==========  ===========  ==========  ===================================

Optional inputs
***************

==========================  =====================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
name                        type                   prefix                     position  documentation
==========================  =====================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
outputBasename              Optional<Filename>     --OUTPUT                             Basename for the three metrics files that are to be written. Resulting files will be:(1) .genotype_concordance_summary_metrics, (2) .genotype_concordance_detail_metrics, (3) .genotype_concordance_contingency_metrics.
argumentsFile               Optional<Array<File>>  --arguments_file                 10  read one or more arguments files and add them to the command line
callSample                  Optional<String>       --CALL_SAMPLE                    10  The name of the call sample within the call VCF. Not required if only one sample exists.
ignoreFilterStatus          Optional<Boolean>      --IGNORE_FILTER_STATUS               Default is false. If true, filter status of sites will be ignored so that we include filtered sites when calculating genotype concordance.
intersectIntervals          Optional<Boolean>      --INTERSECT_INTERVALS                If true, multiple interval lists will be intersected. If false multiple lists will be unioned.
intervals                   Optional<Array<VCF>>   --INTERVALS                          One or more interval list files that will be used to limit the genotype concordance. Note - if intervals are specified, the VCF files must be indexed.
minDP                       Optional<Float>        --MIN_DP                             Genotypes below this depth will have genotypes classified as LowDp.
minGQ                       Optional<Float>        --MIN_GQ                             Genotypes below this genotype quality will have genotypes classified as LowGq.
treatMissingSitesAsHomeRef  Optional<Boolean>      --MISSING_SITES_HOM_REF              Default is false, which follows the GA4GH Scheme. If true, missing sites in the truth
                                                                                        set will be treated as HOM_REF sites and sites missing in both the truth and call sets will be true negatives. Useful when hom ref sites are left out of the truth set. This flag can only be used with a high confidence interval list.
outputAllRows               Optional<Boolean>      --OUTPUT_ALL_ROWS                    If true, output all rows in detailed statistics even when count == 0. When false only output rows with non-zero counts.
outputVcf                   Optional<Boolean>      --OUTPUT_VCF                         Output a VCF annotated with concordance information.
truthSample                 Optional<String>       --TRUTH_SAMPLE                       The name of the truth sample within the truth VCF. Not required if only one sample exists.
useVcfIndex                 Optional<Boolean>      --USE_VCF_INDEX                      If true, use the VCF index, else iterate over the entire VCF
compressionLevel            Optional<Integer>      --COMPRESSION_LEVEL              11  Compression level for all compressed files created (e.g. BAM and GELI).
createIndex                 Optional<Boolean>      --CREATE_INDEX                   11  Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File               Optional<Boolean>      --CREATE_MD5_FILE                11  Whether to create an MD5 digest for any BAM or FASTQ files created.
maxRecordsInRam             Optional<Integer>      --MAX_RECORDS_IN_RAM             11  When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
quiet                       Optional<Boolean>      --QUIET                          11  Whether to suppress job-summary info on System.err.
reference                   Optional<File>         --REFERENCE=SEQUENCE             11  Reference sequence file.
tmpDir                      Optional<String>       --TMP_DIR                        11  Undocumented option
useJdkDeflater              Optional<Boolean>      --use_jdk_deflater               11  Whether to use the JdkDeflater (as opposed to IntelDeflater)
useJdkInflater              Optional<Boolean>      --use_jdk_inflater               11  Whether to use the JdkInflater (as opposed to IntelInflater)
validationStringency        Optional<String>       --VALIDATION_STRINGENCY          11  Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
verbosity                   Optional<String>       --verbosity                      11  The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
==========================  =====================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: Michael Franklin


*GATK4: Genotype Concordance was last updated on 2018-12-24*.
*This page was automatically generated on 2019-04-18*.
