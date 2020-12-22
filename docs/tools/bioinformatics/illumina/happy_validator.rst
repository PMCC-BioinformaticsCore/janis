:orphan:

Hap.py validation
===================================

``happy_validator`` · *1 contributor · 1 version*

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


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.illumina.happy.versions import HapPyValidator_0_3_9

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "happy_validator_step",
           HapPyValidator_0_3_9(
               truthVCF=None,
               compareVCF=None,
               reference=None,
           )
       )
       wf.output("extended", source=happy_validator_step.extended)
       wf.output("summary", source=happy_validator_step.summary)
       wf.output("metrics", source=happy_validator_step.metrics)
       wf.output("vcf", source=happy_validator_step.vcf)
       wf.output("runinfo", source=happy_validator_step.runinfo)
       wf.output("rocOut", source=happy_validator_step.rocOut)
       wf.output("indelLocations", source=happy_validator_step.indelLocations)
       wf.output("indelPassLocations", source=happy_validator_step.indelPassLocations)
       wf.output("snpLocations", source=happy_validator_step.snpLocations)
       wf.output("snpPassLocations", source=happy_validator_step.snpPassLocations)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for happy_validator:

.. code-block:: bash

   # user inputs
   janis inputs happy_validator > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       compareVCF: compareVCF.vcf
       reference: reference.fasta
       truthVCF: truthVCF.vcf




5. Run happy_validator with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       happy_validator





Information
------------

:ID: ``happy_validator``
:URL: *No URL to the documentation was provided*
:Versions: v0.3.9
:Container: pkrusche/hap.py:v0.3.9
:Authors: Michael Franklin
:Citations: None
:Created: 2019-05-15
:Updated: 2019-05-15


Outputs
-----------

==================  ============  ===============
name                type          documentation
==================  ============  ===============
extended            csv
summary             csv
metrics             File
vcf                 Gzipped<VCF>
runinfo             jsonFile
rocOut              File
indelLocations      File
indelPassLocations  File
snpLocations        File
snpPassLocations    File
==================  ============  ===============


Additional configuration (inputs)
---------------------------------

========================  ==================  ============================  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                      type                prefix                          position  documentation
========================  ==================  ============================  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
truthVCF                  VCF                                                        1
compareVCF                VCF                                                        2
reference                 FastaWithIndexes    --reference                               (-r)  Specify a reference file.
reportPrefix              Optional<Filename>  --report-prefix                           (-o)  Filename prefix for report output.
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
engine                    Optional<String>    --engine                                  {xcmp,vcfeval,scmp-somatic,scmp-distance} Comparison engine to use.
engineVcfevalTemplate     Optional<String>    --engine-vcfeval-template                 Vcfeval needs the reference sequence formatted in its own file format (SDF -- run rtg format -o ref.SDF ref.fa). You can specify this here to save time when running hap.py with vcfeval. If no SDF folder is specified, hap.py will create a temporary one.
scmpDistance              Optional<Integer>   --scmp-distance                           For distance-based matching, this is the distance between variants to use.
logfile                   Optional<Filename>  --logfile                                 Write logging information into file rather than to stderr
verbose                   Optional<Boolean>   --verbose                                 Raise logging level from warning to info.
quiet                     Optional<Boolean>   --quiet                                   Set logging level to output errors only.
========================  ==================  ============================  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task happy_validator {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File truthVCF
       File compareVCF
       String? reportPrefix
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       File? intervals
       Boolean? version
       String? scratchPrefix
       String? keepScratch
       File? falsePositives
       File? stratification
       String? stratificationRegion
       String? stratificationFixchr
       Boolean? writeVcf
       Boolean? writeCounts
       Boolean? noWriteCounts
       Boolean? outputVtc
       Boolean? preserveInfo
       String? roc
       Boolean? noRoc
       String? rocRegions
       String? rocFilter
       Int? rocDelta
       Int? ciAlpha
       Boolean? noJson
       Boolean? passOnly
       Boolean? restrictRegions
       Boolean? leftshift
       Boolean? noLeftshift
       Boolean? decompose
       Boolean? noDecompose
       Boolean? bcftoolsNorm
       Boolean? fixchr
       Boolean? noFixchr
       Boolean? bcf
       Boolean? somatic
       Boolean? setGT
       String? gender
       Boolean? preprocessTruth
       Boolean? usefilteredTruth
       Boolean? preprocessingWindowSize
       Boolean? adjustConfRegions
       Boolean? noAdjustConfRegions
       Boolean? noHaplotypeComparison
       Int? windowSize
       Int? xcmpEnumerationThreshold
       String? xcmpExpandHapblocks
       Int? threads
       String? engine
       String? engineVcfevalTemplate
       Int? scmpDistance
       String? logfile
       Boolean? verbose
       Boolean? quiet
     }
     command <<<
       set -e
       /opt/hap.py/bin/hap.py \
         --report-prefix '~{select_first([reportPrefix, "generated"])}' \
         --reference '~{reference}' \
         ~{if defined(intervals) then ("--target-regions '" + intervals + "'") else ""} \
         ~{if (defined(version) && select_first([version])) then "--version" else ""} \
         ~{if defined(scratchPrefix) then ("--scratch-prefix '" + scratchPrefix + "'") else ""} \
         ~{if defined(keepScratch) then ("--keep-scratch '" + keepScratch + "'") else ""} \
         ~{if defined(falsePositives) then ("--false-positives '" + falsePositives + "'") else ""} \
         ~{if defined(stratification) then ("--stratification '" + stratification + "'") else ""} \
         ~{if defined(stratificationRegion) then ("--stratification-region '" + stratificationRegion + "'") else ""} \
         ~{if defined(stratificationFixchr) then ("--stratification-fixchr '" + stratificationFixchr + "'") else ""} \
         ~{if (defined(writeVcf) && select_first([writeVcf])) then "--write-vcf" else ""} \
         ~{if (defined(writeCounts) && select_first([writeCounts])) then "--write-counts" else ""} \
         ~{if (defined(noWriteCounts) && select_first([noWriteCounts])) then "--no-write-counts" else ""} \
         ~{if (defined(outputVtc) && select_first([outputVtc])) then "--output-vtc" else ""} \
         ~{if (defined(preserveInfo) && select_first([preserveInfo])) then "--preserve-info" else ""} \
         ~{if defined(roc) then ("--roc '" + roc + "'") else ""} \
         ~{if (defined(noRoc) && select_first([noRoc])) then "--no-roc" else ""} \
         ~{if defined(rocRegions) then ("--roc-regions '" + rocRegions + "'") else ""} \
         ~{if defined(rocFilter) then ("--roc-filter '" + rocFilter + "'") else ""} \
         ~{if defined(rocDelta) then ("--roc-delta " + rocDelta) else ''} \
         ~{if defined(ciAlpha) then ("--ci-alpha " + ciAlpha) else ''} \
         ~{if (defined(noJson) && select_first([noJson])) then "--no-json" else ""} \
         ~{if (defined(passOnly) && select_first([passOnly])) then "--pass-only" else ""} \
         ~{if (defined(restrictRegions) && select_first([restrictRegions])) then "--restrict-regions" else ""} \
         ~{if (defined(leftshift) && select_first([leftshift])) then "--leftshift" else ""} \
         ~{if (defined(noLeftshift) && select_first([noLeftshift])) then "--no-leftshift" else ""} \
         ~{if (defined(decompose) && select_first([decompose])) then "--decompose" else ""} \
         ~{if (defined(noDecompose) && select_first([noDecompose])) then "--no-decompose" else ""} \
         ~{if (defined(bcftoolsNorm) && select_first([bcftoolsNorm])) then "--bcftools-norm" else ""} \
         ~{if (defined(fixchr) && select_first([fixchr])) then "--fixchr" else ""} \
         ~{if (defined(noFixchr) && select_first([noFixchr])) then "--no-fixchr" else ""} \
         ~{if (defined(bcf) && select_first([bcf])) then "--bcf" else ""} \
         ~{if (defined(somatic) && select_first([somatic])) then "--somatic" else ""} \
         ~{if (defined(setGT) && select_first([setGT])) then "--set-gt" else ""} \
         ~{if defined(gender) then ("--gender '" + gender + "'") else ""} \
         ~{if (defined(preprocessTruth) && select_first([preprocessTruth])) then "--preprocess-truth" else ""} \
         ~{if (defined(usefilteredTruth) && select_first([usefilteredTruth])) then "--usefiltered-truth" else ""} \
         ~{if (defined(preprocessingWindowSize) && select_first([preprocessingWindowSize])) then "--preprocessing-window-size" else ""} \
         ~{if (defined(adjustConfRegions) && select_first([adjustConfRegions])) then "--adjust-conf-regions" else ""} \
         ~{if (defined(noAdjustConfRegions) && select_first([noAdjustConfRegions])) then "--no-adjust-conf-regions" else ""} \
         ~{if (defined(noHaplotypeComparison) && select_first([noHaplotypeComparison])) then "--no-haplotype-comparison" else ""} \
         ~{if defined(windowSize) then ("--window-size " + windowSize) else ''} \
         ~{if defined(xcmpEnumerationThreshold) then ("--xcmp-enumeration-threshold " + xcmpEnumerationThreshold) else ''} \
         ~{if defined(xcmpExpandHapblocks) then ("--xcmp-expand-hapblocks '" + xcmpExpandHapblocks + "'") else ""} \
         ~{if defined(select_first([threads, select_first([runtime_cpu, 1])])) then ("--threads " + select_first([threads, select_first([runtime_cpu, 1])])) else ''} \
         ~{if defined(engine) then ("--engine '" + engine + "'") else ""} \
         ~{if defined(engineVcfevalTemplate) then ("--engine-vcfeval-template '" + engineVcfevalTemplate + "'") else ""} \
         ~{if defined(scmpDistance) then ("--scmp-distance " + scmpDistance) else ''} \
         --logfile '~{select_first([logfile, "generated--log.txt"])}' \
         ~{if (defined(verbose) && select_first([verbose])) then "--verbose" else ""} \
         ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
         '~{truthVCF}' \
         '~{compareVCF}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 2, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "pkrusche/hap.py:v0.3.9"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File extended = (select_first([reportPrefix, "generated"]) + ".extended.csv")
       File summary = (select_first([reportPrefix, "generated"]) + ".summary.csv")
       File metrics = (select_first([reportPrefix, "generated"]) + ".metrics.json.gz")
       File vcf = (select_first([reportPrefix, "generated"]) + ".vcf.gz")
       File vcf_tbi = (select_first([reportPrefix, "generated"]) + ".vcf.gz") + ".tbi"
       File runinfo = (select_first([reportPrefix, "generated"]) + ".runinfo.json")
       File rocOut = (select_first([reportPrefix, "generated"]) + ".roc.all.csv.gz")
       File indelLocations = (select_first([reportPrefix, "generated"]) + ".roc.Locations.INDEL.csv.gz")
       File indelPassLocations = (select_first([reportPrefix, "generated"]) + ".roc.Locations.INDEL.PASS.csv.gz")
       File snpLocations = (select_first([reportPrefix, "generated"]) + ".roc.Locations.SNP.csv.gz")
       File snpPassLocations = (select_first([reportPrefix, "generated"]) + ".roc.Locations.SNP.PASS.csv.gz")
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Hap.py validation
   doc: |-
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

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: pkrusche/hap.py:v0.3.9

   inputs:
   - id: truthVCF
     label: truthVCF
     type: File
     inputBinding:
       position: 1
   - id: compareVCF
     label: compareVCF
     type: File
     inputBinding:
       position: 2
   - id: reportPrefix
     label: reportPrefix
     doc: (-o)  Filename prefix for report output.
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --report-prefix
   - id: reference
     label: reference
     doc: (-r)  Specify a reference file.
     type: File
     secondaryFiles:
     - pattern: .fai
     - pattern: .amb
     - pattern: .ann
     - pattern: .bwt
     - pattern: .pac
     - pattern: .sa
     - pattern: ^.dict
     inputBinding:
       prefix: --reference
   - id: intervals
     label: intervals
     doc: (-T)  Restrict analysis to given (dense) regions (using -T in bcftools).
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --target-regions
   - id: version
     label: version
     doc: (-v) Show version number and exit.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --version
   - id: scratchPrefix
     label: scratchPrefix
     doc: Directory for scratch files.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --scratch-prefix
   - id: keepScratch
     label: keepScratch
     doc: Filename prefix for scratch report output. Annotation format in input VCF file.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --keep-scratch
   - id: falsePositives
     label: falsePositives
     doc: |-
       (-f)  False positive / confident call regions (.bed or .bed.gz). Calls outside these regions will be labelled as UNK.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --false-positives
   - id: stratification
     label: stratification
     doc: |2-
        Stratification file list (TSV format -- first column is region name, second column is file name).
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --stratification
   - id: stratificationRegion
     label: stratificationRegion
     doc: Add single stratification region, e.g. --stratification-region TEST:test.bed
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --stratification-region
   - id: stratificationFixchr
     label: stratificationFixchr
     doc: ' Add chr prefix to stratification files if necessary'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --stratification-fixchr
   - id: writeVcf
     label: writeVcf
     doc: (-V) Write an annotated VCF.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --write-vcf
   - id: writeCounts
     label: writeCounts
     doc: (-X) Write advanced counts and metrics.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --write-counts
   - id: noWriteCounts
     label: noWriteCounts
     doc: Do not write advanced counts and metrics.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-write-counts
   - id: outputVtc
     label: outputVtc
     doc: |-
       Write VTC field in the final VCF which gives the counts each position has contributed to.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --output-vtc
   - id: preserveInfo
     label: preserveInfo
     doc: |-
       When using XCMP, preserve and merge the INFO fields in truth and query. Useful for ROC computation.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --preserve-info
   - id: roc
     label: roc
     doc: Select a feature to produce a ROC on (INFO feature, QUAL, GQX, ...).
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --roc
   - id: noRoc
     label: noRoc
     doc: |-
       Disable ROC computation and only output summary statistics for more concise output.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-roc
   - id: rocRegions
     label: rocRegions
     doc: |2-
        Select a list of regions to compute ROCs in. By default, only the '*' region will produce ROC output (aggregate variant counts).
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --roc-regions
   - id: rocFilter
     label: rocFilter
     doc: ' Select a filter to ignore when making ROCs.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --roc-filter
   - id: rocDelta
     label: rocDelta
     doc: ' Minimum spacing between ROC QQ levels.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --roc-delta
   - id: ciAlpha
     label: ciAlpha
     doc: |-
       Confidence level for Jeffrey's CI for recall, precision and fraction of non-assessed calls.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --ci-alpha
   - id: noJson
     label: noJson
     doc: Disable JSON file output.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-json
   - id: passOnly
     label: passOnly
     doc: Keep only PASS variants.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --pass-only
   - id: restrictRegions
     label: restrictRegions
     doc: (-R)  Restrict analysis to given (sparse) regions (using -R in bcftools).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --restrict-regions
   - id: leftshift
     label: leftshift
     doc: (-L) Left-shift variants safely.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --leftshift
   - id: noLeftshift
     label: noLeftshift
     doc: Do not left-shift variants safely.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-leftshift
   - id: decompose
     label: decompose
     doc: Decompose variants into primitives. This results in more granular counts.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --decompose
   - id: noDecompose
     label: noDecompose
     doc: (-D) Do not decompose variants into primitives.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-decompose
   - id: bcftoolsNorm
     label: bcftoolsNorm
     doc: |-
       Enable preprocessing through bcftools norm -c x -D (requires external preprocessing to be switched on).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --bcftools-norm
   - id: fixchr
     label: fixchr
     doc: |-
       Add chr prefix to VCF records where necessary (default: auto, attempt to match reference).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --fixchr
   - id: noFixchr
     label: noFixchr
     doc: |-
       Do not add chr prefix to VCF records (default: auto, attempt to match reference).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-fixchr
   - id: bcf
     label: bcf
     doc: |-
       Use BCF internally. This is the default when the input file is in BCF format already. Using BCF can speed up temp file access, but may fail for VCF files that have broken headers or records that don't comply with the header.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --bcf
   - id: somatic
     label: somatic
     doc: |-
       Assume the input file is a somatic call file and squash all columns into one, putting all FORMATs into INFO + use half genotypes (see also --set-gt). This will replace all sample columns and replace them with a single one. This is used to treat Strelka somatic files Possible values for this parameter: half / hemi / het / hom / half to assign one of the following genotypes to the resulting sample: 1 | 0/1 | 1/1 | ./1. This will replace all sample columns and replace them with a single one.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --somatic
   - id: setGT
     label: setGT
     doc: |-
       This is used to treat Strelka somatic files Possible values for this parameter: half / hemi / het / hom / half to assign one of the following genotypes to the resulting sample: 1 | 0/1 | 1/1 | ./1. This will replace all sample columns and replace them with a single one.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --set-gt
   - id: gender
     label: gender
     doc: |-
       ({male,female,auto,none})  Specify gender. This determines how haploid calls on chrX get treated: for male samples, all non-ref calls (in the truthset only when running through hap.py) are given a 1/1 genotype.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --gender
   - id: preprocessTruth
     label: preprocessTruth
     doc: |-
       Preprocess truth file with same settings as query (default is to accept truth in original format).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --preprocess-truth
   - id: usefilteredTruth
     label: usefilteredTruth
     doc: |-
       Use filtered variant calls in truth file (by default, only PASS calls in the truth file are used)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --usefiltered-truth
   - id: preprocessingWindowSize
     label: preprocessingWindowSize
     doc: |2-
        Preprocessing window size (variants further apart than that size are not expected to interfere).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --preprocessing-window-size
   - id: adjustConfRegions
     label: adjustConfRegions
     doc: |2-
        Adjust confident regions to include variant locations. Note this will only include variants that are included in the CONF regions already when viewing with bcftools; this option only makes sure insertions are padded correctly in the CONF regions (to capture these, both the base before and after must be contained in the bed file).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --adjust-conf-regions
   - id: noAdjustConfRegions
     label: noAdjustConfRegions
     doc: ' Do not adjust confident regions for insertions.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-adjust-conf-regions
   - id: noHaplotypeComparison
     label: noHaplotypeComparison
     doc: (--unhappy)  Disable haplotype comparison (only count direct GT matches as
       TP).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-haplotype-comparison
   - id: windowSize
     label: windowSize
     doc: |-
       (-w)  Minimum distance between variants such that they fall into the same superlocus.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --window-size
   - id: xcmpEnumerationThreshold
     label: xcmpEnumerationThreshold
     doc: ' Enumeration threshold / maximum number of sequences to enumerate per block.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --xcmp-enumeration-threshold
   - id: xcmpExpandHapblocks
     label: xcmpExpandHapblocks
     doc: ' Expand haplotype blocks by this many basepairs left and right.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --xcmp-expand-hapblocks
   - id: threads
     label: threads
     doc: Number of threads to use. Comparison engine to use.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --threads
       valueFrom: |-
         $([inputs.runtime_cpu, 2, 1].filter(function (inner) { return inner != null })[0])
   - id: engine
     label: engine
     doc: ' {xcmp,vcfeval,scmp-somatic,scmp-distance} Comparison engine to use.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --engine
   - id: engineVcfevalTemplate
     label: engineVcfevalTemplate
     doc: |2-
        Vcfeval needs the reference sequence formatted in its own file format (SDF -- run rtg format -o ref.SDF ref.fa). You can specify this here to save time when running hap.py with vcfeval. If no SDF folder is specified, hap.py will create a temporary one.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --engine-vcfeval-template
   - id: scmpDistance
     label: scmpDistance
     doc: ' For distance-based matching, this is the distance between variants to use.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --scmp-distance
   - id: logfile
     label: logfile
     doc: Write logging information into file rather than to stderr
     type:
     - string
     - 'null'
     default: generated--log.txt
     inputBinding:
       prefix: --logfile
   - id: verbose
     label: verbose
     doc: Raise logging level from warning to info.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --verbose
   - id: quiet
     label: quiet
     doc: Set logging level to output errors only.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --quiet

   outputs:
   - id: extended
     label: extended
     type: File
     outputBinding:
       glob: $((inputs.reportPrefix + ".extended.csv"))
       outputEval: $((inputs.reportPrefix.basename + ".extended.csv"))
       loadContents: false
   - id: summary
     label: summary
     type: File
     outputBinding:
       glob: $((inputs.reportPrefix + ".summary.csv"))
       outputEval: $((inputs.reportPrefix.basename + ".summary.csv"))
       loadContents: false
   - id: metrics
     label: metrics
     type: File
     outputBinding:
       glob: $((inputs.reportPrefix + ".metrics.json.gz"))
       outputEval: $((inputs.reportPrefix.basename + ".metrics.json.gz"))
       loadContents: false
   - id: vcf
     label: vcf
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $((inputs.reportPrefix + ".vcf.gz"))
       outputEval: $((inputs.reportPrefix.basename + ".vcf.gz"))
       loadContents: false
   - id: runinfo
     label: runinfo
     type: File
     outputBinding:
       glob: $((inputs.reportPrefix + ".runinfo.json"))
       outputEval: $((inputs.reportPrefix.basename + ".runinfo.json"))
       loadContents: false
   - id: rocOut
     label: rocOut
     type: File
     outputBinding:
       glob: $((inputs.reportPrefix + ".roc.all.csv.gz"))
       outputEval: $((inputs.reportPrefix.basename + ".roc.all.csv.gz"))
       loadContents: false
   - id: indelLocations
     label: indelLocations
     type: File
     outputBinding:
       glob: $((inputs.reportPrefix + ".roc.Locations.INDEL.csv.gz"))
       outputEval: $((inputs.reportPrefix.basename + ".roc.Locations.INDEL.csv.gz"))
       loadContents: false
   - id: indelPassLocations
     label: indelPassLocations
     type: File
     outputBinding:
       glob: $((inputs.reportPrefix + ".roc.Locations.INDEL.PASS.csv.gz"))
       outputEval: $((inputs.reportPrefix.basename + ".roc.Locations.INDEL.PASS.csv.gz"))
       loadContents: false
   - id: snpLocations
     label: snpLocations
     type: File
     outputBinding:
       glob: $((inputs.reportPrefix + ".roc.Locations.SNP.csv.gz"))
       outputEval: $((inputs.reportPrefix.basename + ".roc.Locations.SNP.csv.gz"))
       loadContents: false
   - id: snpPassLocations
     label: snpPassLocations
     type: File
     outputBinding:
       glob: $((inputs.reportPrefix + ".roc.Locations.SNP.PASS.csv.gz"))
       outputEval: $((inputs.reportPrefix.basename + ".roc.Locations.SNP.PASS.csv.gz"))
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: /opt/hap.py/bin/hap.py
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: happy_validator


