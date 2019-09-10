
.. include:: happy_validator_v0.3.9

Hap.py validation
===================================

Description
-------------

Tool identifier: ``happy_validator``

Tool path: ``janis_bioinformatics.tools.illumina.happy.versions import HapPyValidator_0_3_9``

Version: v0.3.9

Container: ``pkrusche/hap.py:v0.3.9``



Documentation
-------------

URL
******
*No URL to the documentation was provided*

Tool documentation
******************
usage: Haplotype Comparison 
    [-h] [-v] [-r REF] [-o REPORTS_PREFIX]
    [--scratch-prefix SCRATCH_PREFIX] [--keep-scratch]
    [-t {xcmp,ga4gh}] [-f FP_BEDFILE]
    [--stratification STRAT_TSV]
    [--stratification-region STRAT_REGIONS]
    [--stratification-fixchr] [-V] [-X]
    [--no-write-counts] [--output-vtc]
    [--preserve-info] [--roc ROC] [--no-roc]
    [--roc-regions ROC_REGIONS]
    [--roc-filter ROC_FILTER] [--roc-delta ROC_DELTA]
    [--ci-alpha CI_ALPHA] [--no-json]
    [--location LOCATIONS] [--pass-only]
    [--filters-only FILTERS_ONLY] [-R REGIONS_BEDFILE]
    [-T TARGETS_BEDFILE] [-L] [--no-leftshift]
    [--decompose] [-D] [--bcftools-norm] [--fixchr]
    [--no-fixchr] [--bcf] [--somatic]
    [--set-gt {half,hemi,het,hom,first}]
    [--gender {male,female,auto,none}]
    [--preprocess-truth] [--usefiltered-truth]
    [--preprocessing-window-size PREPROCESS_WINDOW]
    [--adjust-conf-regions] [--no-adjust-conf-regions]
    [--unhappy] [-w WINDOW]
    [--xcmp-enumeration-threshold MAX_ENUM]
    [--xcmp-expand-hapblocks HB_EXPAND]
    [--threads THREADS]
    [--engine {xcmp,vcfeval,scmp-somatic,scmp-distance}]
    [--engine-vcfeval-path ENGINE_VCFEVAL]
    [--engine-vcfeval-template ENGINE_VCFEVAL_TEMPLATE]
    [--scmp-distance ENGINE_SCMP_DISTANCE]
    [--logfile LOGFILE] [--verbose | --quiet]
    [_vcfs [_vcfs ...]]
positional arguments:
  _vcfs                 Two VCF files.

Outputs
-------
==================  ====================  ===============
name                type                  documentation
==================  ====================  ===============
extended            csv
summary             csv
metrics             File
vcf                 CompressedIndexedVCF
runinfo             jsonFile
rocOut              File
indelLocations      File
indelPassLocations  File
snpLocations        File
snpPassLocations    File
==================  ====================  ===============

Inputs
------
Find the inputs below

Required inputs
***************

============  =============  ===============  ==========  ========================================
name          type           prefix             position  documentation
============  =============  ===============  ==========  ========================================
truthVCF      VCF                                      1
compareVCF    VCF                                      2
reportPrefix  String         --report-prefix              (-o)  Filename prefix for report output.
reference     FastaWithDict  --reference                  (-r)  Specify a reference file.
============  =============  ===============  ==========  ========================================

Optional inputs
***************

========================  ==================  ============================  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                      type                prefix                        position    documentation
========================  ==================  ============================  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
intervals                 Optional<bed>       --target-regions                          (-T)  Restrict analysis to given (dense) regions (using -T in bcftools).
version                   Optional<Boolean>   --version                                 (-v) Show version number and exit.
scratchPrefix             Optional<String>    --scratch-prefix                          Directory for scratch files.
keepScratch               Optional<String>    --keep-scratch                            Filename prefix for scratch report output. Annotation format in input VCF file.
falsePositives            Optional<bed>       --false-positives                         (-f)  False positive / confident call regions (.bed or .bed.gz). Calls outside these regions will be labelled as UNK.
stratification            Optional<tsv>       --stratification                          Stratification file list (TSV format -- first column is region name, second column is file name).
stratificationRegion      Optional<String>    --stratification-region                   Add single stratification region, e.g. --stratification-region TEST:test.bed
stratificationFixchr      Optional<String>    --stratification-fixchr                   Add chr prefix to stratification files if necessary
writeVcf                  Optional<Boolean>   --write-vcf                               (-V) Write an annotated VCF.
writeCounts               Optional<Boolean>   --write-counts                            (-X) Write advanced counts and metrics.
noWriteCounts             Optional<Boolean>   --no-write-counts                         Do not write advanced counts and metrics.
outputVtc                 Optional<Boolean>   --output-vtc                              Write VTC field in the final VCF which gives the counts each position has contributed to.
preserveInfo              Optional<Boolean>   --preserve-info                           When using XCMP, preserve and merge the INFO fields in truth and query. Useful for ROC computation.
roc                       Optional<String>    --roc                                     Select a feature to produce a ROC on (INFO feature, QUAL, GQX, ...).
noRoc                     Optional<Boolean>   --no-roc                                  Disable ROC computation and only output summary statistics for more concise output.
rocRegions                Optional<String>    --roc-regions                             Select a list of regions to compute ROCs in. By default, only the '*' region will produce ROC output (aggregate variant counts).
rocFilter                 Optional<String>    --roc-filter                              Select a filter to ignore when making ROCs.
rocDelta                  Optional<Integer>   --roc-delta                               Minimum spacing between ROC QQ levels.
ciAlpha                   Optional<Integer>   --ci-alpha                                Confidence level for Jeffrey's CI for recall, precision and fraction of non-assessed calls.
noJson                    Optional<Boolean>   --no-json                                 Disable JSON file output.
passOnly                  Optional<Boolean>   --pass-only                               Keep only PASS variants.
restrictRegions           Optional<Boolean>   --restrict-regions                        (-R)  Restrict analysis to given (sparse) regions (using -R in bcftools).
leftshift                 Optional<Boolean>   --leftshift                               (-L) Left-shift variants safely.
noLeftshift               Optional<Boolean>   --no-leftshift                            Do not left-shift variants safely.
decompose                 Optional<Boolean>   --decompose                               Decompose variants into primitives. This results in more granular counts.
noDecompose               Optional<Boolean>   --no-decompose                            (-D) Do not decompose variants into primitives.
bcftoolsNorm              Optional<Boolean>   --bcftools-norm                           Enable preprocessing through bcftools norm -c x -D (requires external preprocessing to be switched on).
fixchr                    Optional<Boolean>   --fixchr                                  Add chr prefix to VCF records where necessary (default: auto, attempt to match reference).
noFixchr                  Optional<Boolean>   --no-fixchr                               Do not add chr prefix to VCF records (default: auto, attempt to match reference).
bcf                       Optional<Boolean>   --bcf                                     Use BCF internally. This is the default when the input file is in BCF format already. Using BCF can speed up temp file access, but may fail for VCF files that have broken headers or records that don't comply with the header.
somatic                   Optional<Boolean>   --somatic                                 Assume the input file is a somatic call file and squash all columns into one, putting all FORMATs into INFO + use half genotypes (see also --set-gt). This will replace all sample columns and replace them with a single one. This is used to treat Strelka somatic files Possible values for this parameter: half / hemi / het / hom / half to assign one of the following genotypes to the resulting sample: 1 | 0/1 | 1/1 | ./1. This will replace all sample columns and replace them with a single one.
setGT                     Optional<Boolean>   --set-gt                                  This is used to treat Strelka somatic files Possible values for this parameter: half / hemi / het / hom / half to assign one of the following genotypes to the resulting sample: 1 | 0/1 | 1/1 | ./1. This will replace all sample columns and replace them with a single one.
gender                    Optional<String>    --gender                                  ({male,female,auto,none})  Specify gender. This determines how haploid calls on chrX get treated: for male samples, all non-ref calls (in the truthset only when running through hap.py) are given a 1/1 genotype.
preprocessTruth           Optional<Boolean>   --preprocess-truth                        Preprocess truth file with same settings as query (default is to accept truth in original format).
usefilteredTruth          Optional<Boolean>   --usefiltered-truth                       Use filtered variant calls in truth file (by default, only PASS calls in the truth file are used)
preprocessingWindowSize   Optional<Boolean>   --preprocessing-window-size               Preprocessing window size (variants further apart than that size are not expected to interfere).
adjustConfRegions         Optional<Boolean>   --adjust-conf-regions                     Adjust confident regions to include variant locations. Note this will only include variants that are included in the CONF regions already when viewing with bcftools; this option only makes sure insertions are padded correctly in the CONF regions (to capture these, both the base before and after must be contained in the bed file).
noAdjustConfRegions       Optional<Boolean>   --no-adjust-conf-regions                  Do not adjust confident regions for insertions.
noHaplotypeComparison     Optional<Boolean>   --no-haplotype-comparison                 (--unhappy)  Disable haplotype comparison (only count direct GT matches as TP).
windowSize                Optional<Integer>   --window-size                             (-w)  Minimum distance between variants such that they fall into the same superlocus.
xcmpEnumerationThreshold  Optional<Integer>   --xcmp-enumeration-threshold              Enumeration threshold / maximum number of sequences to enumerate per block.
xcmpExpandHapblocks       Optional<String>    --xcmp-expand-hapblocks                   Expand haplotype blocks by this many basepairs left and right.
threads                   Optional<Integer>   --threads                                 Number of threads to use. Comparison engine to use.
engineVcfevalTemplate     Optional<String>    --engine-vcfeval-template                 Vcfeval needs the reference sequence formatted in its own file format (SDF -- run rtg format -o ref.SDF ref.fa). You can specify this here to save time when running hap.py with vcfeval. If no SDF folder is specified, hap.py will create a temporary one.
scmpDistance              Optional<Integer>   --scmp-distance                           For distance-based matching, this is the distance between variants to use.
logfile                   Optional<Filename>  --logfile                                 Write logging information into file rather than to stderr
verbose                   Optional<Boolean>   --verbose                                 Raise logging level from warning to info.
quiet                     Optional<Boolean>   --quiet                                   Set logging level to output errors only.
========================  ==================  ============================  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: **Unknown**


*Hap.py validation was last updated on 2019-05-15*.
*This page was automatically generated on 2019-09-10*.
