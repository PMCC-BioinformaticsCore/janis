:orphan:

Hap.py validation
===================================

1 contributor · 1 version

:ID: ``happy_validator``
:Python: ``janis_bioinformatics.tools.illumina.happy.versions import HapPyValidator_0_3_9``
:Versions: v0.3.9
:Container: pkrusche/hap.py:v0.3.9
:Authors: Michael Franklin
:Citations: None
:Created: 2019-05-15
:Updated: 2019-05-15
:Required inputs:
   - ``truthVCF: VCF``

   - ``compareVCF: VCF``

   - ``reference: FastaWithDict``
:Outputs: 
   - ``extended: csv``

   - ``summary: csv``

   - ``metrics: File``

   - ``vcf: CompressedIndexedVCF``

   - ``runinfo: jsonFile``

   - ``rocOut: File``

   - ``indelLocations: File``

   - ``indelPassLocations: File``

   - ``snpLocations: File``

   - ``snpPassLocations: File``

Documentation
-------------

URL: *No URL to the documentation was provided*

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

------

Additional configuration (inputs)
---------------------------------

========================  ==================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                      type                documentation
========================  ==================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
truthVCF                  VCF
compareVCF                VCF
reference                 FastaWithDict       (-r)  Specify a reference file.
reportPrefix              Optional<Filename>  (-o)  Filename prefix for report output.
intervals                 Optional<bed>       (-T)  Restrict analysis to given (dense) regions (using -T in bcftools).
version                   Optional<Boolean>   (-v) Show version number and exit.
scratchPrefix             Optional<String>    Directory for scratch files.
keepScratch               Optional<String>    Filename prefix for scratch report output. Annotation format in input VCF file.
falsePositives            Optional<bed>       (-f)  False positive / confident call regions (.bed or .bed.gz). Calls outside these regions will be labelled as UNK.
stratification            Optional<tsv>       Stratification file list (TSV format -- first column is region name, second column is file name).
stratificationRegion      Optional<String>    Add single stratification region, e.g. --stratification-region TEST:test.bed
stratificationFixchr      Optional<String>    Add chr prefix to stratification files if necessary
writeVcf                  Optional<Boolean>   (-V) Write an annotated VCF.
writeCounts               Optional<Boolean>   (-X) Write advanced counts and metrics.
noWriteCounts             Optional<Boolean>   Do not write advanced counts and metrics.
outputVtc                 Optional<Boolean>   Write VTC field in the final VCF which gives the counts each position has contributed to.
preserveInfo              Optional<Boolean>   When using XCMP, preserve and merge the INFO fields in truth and query. Useful for ROC computation.
roc                       Optional<String>    Select a feature to produce a ROC on (INFO feature, QUAL, GQX, ...).
noRoc                     Optional<Boolean>   Disable ROC computation and only output summary statistics for more concise output.
rocRegions                Optional<String>    Select a list of regions to compute ROCs in. By default, only the '*' region will produce ROC output (aggregate variant counts).
rocFilter                 Optional<String>    Select a filter to ignore when making ROCs.
rocDelta                  Optional<Integer>   Minimum spacing between ROC QQ levels.
ciAlpha                   Optional<Integer>   Confidence level for Jeffrey's CI for recall, precision and fraction of non-assessed calls.
noJson                    Optional<Boolean>   Disable JSON file output.
passOnly                  Optional<Boolean>   Keep only PASS variants.
restrictRegions           Optional<Boolean>   (-R)  Restrict analysis to given (sparse) regions (using -R in bcftools).
leftshift                 Optional<Boolean>   (-L) Left-shift variants safely.
noLeftshift               Optional<Boolean>   Do not left-shift variants safely.
decompose                 Optional<Boolean>   Decompose variants into primitives. This results in more granular counts.
noDecompose               Optional<Boolean>   (-D) Do not decompose variants into primitives.
bcftoolsNorm              Optional<Boolean>   Enable preprocessing through bcftools norm -c x -D (requires external preprocessing to be switched on).
fixchr                    Optional<Boolean>   Add chr prefix to VCF records where necessary (default: auto, attempt to match reference).
noFixchr                  Optional<Boolean>   Do not add chr prefix to VCF records (default: auto, attempt to match reference).
bcf                       Optional<Boolean>   Use BCF internally. This is the default when the input file is in BCF format already. Using BCF can speed up temp file access, but may fail for VCF files that have broken headers or records that don't comply with the header.
somatic                   Optional<Boolean>   Assume the input file is a somatic call file and squash all columns into one, putting all FORMATs into INFO + use half genotypes (see also --set-gt). This will replace all sample columns and replace them with a single one. This is used to treat Strelka somatic files Possible values for this parameter: half / hemi / het / hom / half to assign one of the following genotypes to the resulting sample: 1 | 0/1 | 1/1 | ./1. This will replace all sample columns and replace them with a single one.
setGT                     Optional<Boolean>   This is used to treat Strelka somatic files Possible values for this parameter: half / hemi / het / hom / half to assign one of the following genotypes to the resulting sample: 1 | 0/1 | 1/1 | ./1. This will replace all sample columns and replace them with a single one.
gender                    Optional<String>    ({male,female,auto,none})  Specify gender. This determines how haploid calls on chrX get treated: for male samples, all non-ref calls (in the truthset only when running through hap.py) are given a 1/1 genotype.
preprocessTruth           Optional<Boolean>   Preprocess truth file with same settings as query (default is to accept truth in original format).
usefilteredTruth          Optional<Boolean>   Use filtered variant calls in truth file (by default, only PASS calls in the truth file are used)
preprocessingWindowSize   Optional<Boolean>   Preprocessing window size (variants further apart than that size are not expected to interfere).
adjustConfRegions         Optional<Boolean>   Adjust confident regions to include variant locations. Note this will only include variants that are included in the CONF regions already when viewing with bcftools; this option only makes sure insertions are padded correctly in the CONF regions (to capture these, both the base before and after must be contained in the bed file).
noAdjustConfRegions       Optional<Boolean>   Do not adjust confident regions for insertions.
noHaplotypeComparison     Optional<Boolean>   (--unhappy)  Disable haplotype comparison (only count direct GT matches as TP).
windowSize                Optional<Integer>   (-w)  Minimum distance between variants such that they fall into the same superlocus.
xcmpEnumerationThreshold  Optional<Integer>   Enumeration threshold / maximum number of sequences to enumerate per block.
xcmpExpandHapblocks       Optional<String>    Expand haplotype blocks by this many basepairs left and right.
threads                   Optional<Integer>   Number of threads to use. Comparison engine to use.
engine                    Optional<String>    {xcmp,vcfeval,scmp-somatic,scmp-distance} Comparison engine to use.
engineVcfevalTemplate     Optional<String>    Vcfeval needs the reference sequence formatted in its own file format (SDF -- run rtg format -o ref.SDF ref.fa). You can specify this here to save time when running hap.py with vcfeval. If no SDF folder is specified, hap.py will create a temporary one.
scmpDistance              Optional<Integer>   For distance-based matching, this is the distance between variants to use.
logfile                   Optional<Filename>  Write logging information into file rather than to stderr
verbose                   Optional<Boolean>   Raise logging level from warning to info.
quiet                     Optional<Boolean>   Set logging level to output errors only.
========================  ==================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

