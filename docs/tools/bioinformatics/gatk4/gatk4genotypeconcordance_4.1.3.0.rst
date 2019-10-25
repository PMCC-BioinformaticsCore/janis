:orphan:

GATK4: Genotype Concordance
======================================================

1 contributor · 4 versions

:ID: ``Gatk4GenotypeConcordance``
:Python: ``janis_bioinformatics.tools.gatk4.genotypeconcordance.versions import Gatk4GenotypeConcordance_4_1_3``
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.3.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24
:Required inputs:
   - ``callVCF: CompressedIndexedVCF``

   - ``truthVCF: IndexedVCF``
:Outputs: 
   - ``summaryMetrics: File``

   - ``detailMetrics: File``

   - ``contingencyMetrics: File``

Documentation
-------------

URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_vcf_GenotypeConcordance.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_vcf_GenotypeConcordance.php>`_

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

------

Additional configuration (inputs)
---------------------------------

==========================  =====================  ================================================================================================================================================================================================================================================================================================================================================================================================
name                        type                   documentation
==========================  =====================  ================================================================================================================================================================================================================================================================================================================================================================================================
callVCF                     CompressedIndexedVCF   The VCF containing the call sample
truthVCF                    IndexedVCF             The VCF containing the truth sample
outputBasename              Optional<Filename>     Basename for the three metrics files that are to be written. Resulting files will be:(1) .genotype_concordance_summary_metrics, (2) .genotype_concordance_detail_metrics, (3) .genotype_concordance_contingency_metrics.
argumentsFile               Optional<Array<File>>  read one or more arguments files and add them to the command line
callSample                  Optional<String>       The name of the call sample within the call VCF. Not required if only one sample exists.
ignoreFilterStatus          Optional<Boolean>      Default is false. If true, filter status of sites will be ignored so that we include filtered sites when calculating genotype concordance.
intersectIntervals          Optional<Boolean>      If true, multiple interval lists will be intersected. If false multiple lists will be unioned.
intervals                   Optional<Array<VCF>>   One or more interval list files that will be used to limit the genotype concordance. Note - if intervals are specified, the VCF files must be indexed.
minDP                       Optional<Float>        Genotypes below this depth will have genotypes classified as LowDp.
minGQ                       Optional<Float>        Genotypes below this genotype quality will have genotypes classified as LowGq.
treatMissingSitesAsHomeRef  Optional<Boolean>      Default is false, which follows the GA4GH Scheme. If true, missing sites in the truth
                                                   set will be treated as HOM_REF sites and sites missing in both the truth and call sets will be true negatives. Useful when hom ref sites are left out of the truth set. This flag can only be used with a high confidence interval list.
outputAllRows               Optional<Boolean>      If true, output all rows in detailed statistics even when count == 0. When false only output rows with non-zero counts.
outputVcf                   Optional<Boolean>      Output a VCF annotated with concordance information.
truthSample                 Optional<String>       The name of the truth sample within the truth VCF. Not required if only one sample exists.
useVcfIndex                 Optional<Boolean>      If true, use the VCF index, else iterate over the entire VCF
compressionLevel            Optional<Integer>      Compression level for all compressed files created (e.g. BAM and GELI).
createIndex                 Optional<Boolean>      Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File               Optional<Boolean>      Whether to create an MD5 digest for any BAM or FASTQ files created.
maxRecordsInRam             Optional<Integer>      When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
quiet                       Optional<Boolean>      Whether to suppress job-summary info on System.err.
reference                   Optional<File>         Reference sequence file.
tmpDir                      Optional<String>       Undocumented option
useJdkDeflater              Optional<Boolean>      Whether to use the JdkDeflater (as opposed to IntelDeflater)
useJdkInflater              Optional<Boolean>      Whether to use the JdkInflater (as opposed to IntelInflater)
validationStringency        Optional<String>       Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
verbosity                   Optional<String>       The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
==========================  =====================  ================================================================================================================================================================================================================================================================================================================================================================================================

