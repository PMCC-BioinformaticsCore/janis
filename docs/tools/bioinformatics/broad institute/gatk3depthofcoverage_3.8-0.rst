:orphan:

GATK3 DepthOfCoverage: Determine coverage at different levels of partitioning and aggregation.
=====================================================================================================================

``Gatk3DepthOfCoverage`` · *1 contributor · 2 versions*

Overview
This tool processes a set of bam files to determine coverage at different levels of partitioning and aggregation. Coverage can be analyzed per locus, per interval, per gene, or in total; can be partitioned by sample, by read group, by technology, by center, or by library; and can be summarized by mean, median, quartiles, and/or percentage of bases covered to or beyond a threshold. Additionally, reads and bases can be filtered by mapping or base quality score.

Input
One or more bam files (with proper headers) to be analyzed for coverage statistics
(Optional) A REFSEQ file to aggregate coverage to the gene level (for information about creating the REFSEQ Rod, please consult the online documentation)
Output
Tables pertaining to different coverage summaries. Suffix on the table files declares the contents:

no suffix: per locus coverage
_summary: total, mean, median, quartiles, and threshold proportions, aggregated over all bases
_statistics: coverage histograms (# locus with X coverage), aggregated over all bases
_interval_summary: total, mean, median, quartiles, and threshold proportions, aggregated per interval
_interval_statistics: 2x2 table of # of intervals covered to >= X depth in >=Y samples
_gene_summary: total, mean, median, quartiles, and threshold proportions, aggregated per gene
_gene_statistics: 2x2 table of # of genes covered to >= X depth in >= Y samples
_cumulative_coverage_counts: coverage histograms (# locus with >= X coverage), aggregated over all bases
_cumulative_coverage_proportions: proprotions of loci with >= X coverage, aggregated over all bases


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk3.depthofcoverage.versions import GATK3DepthOfCoverage_3_8_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk3depthofcoverage_step",
           GATK3DepthOfCoverage_3_8_0(
               bam=None,
               reference=None,
               outputPrefix=None,
           )
       )
       wf.output("sample", source=gatk3depthofcoverage_step.sample)
       wf.output("sampleCumulativeCoverageCounts", source=gatk3depthofcoverage_step.sampleCumulativeCoverageCounts)
       wf.output("sampleCumulativeCoverageProportions", source=gatk3depthofcoverage_step.sampleCumulativeCoverageProportions)
       wf.output("sampleIntervalStatistics", source=gatk3depthofcoverage_step.sampleIntervalStatistics)
       wf.output("sampleIntervalSummary", source=gatk3depthofcoverage_step.sampleIntervalSummary)
       wf.output("sampleStatistics", source=gatk3depthofcoverage_step.sampleStatistics)
       wf.output("sampleSummary", source=gatk3depthofcoverage_step.sampleSummary)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk3DepthOfCoverage:

.. code-block:: bash

   # user inputs
   janis inputs Gatk3DepthOfCoverage > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       outputPrefix: <value>
       reference: reference.fasta




5. Run Gatk3DepthOfCoverage with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk3DepthOfCoverage





Information
------------

:ID: ``Gatk3DepthOfCoverage``
:URL: `https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tooldocs/3.8-0/org_broadinstitute_gatk_engine_CommandLineGATK.html <https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tooldocs/3.8-0/org_broadinstitute_gatk_engine_CommandLineGATK.html>`_
:Versions: 3.8-1, 3.8-0
:Container: broadinstitute/gatk3:3.8-0
:Authors: Jiaan Yu
:Citations: 
:Created: 2020-04-09
:Updated: 2020-04-09


Outputs
-----------

===================================  ========  ===============
name                                 type      documentation
===================================  ========  ===============
sample                               TextFile
sampleCumulativeCoverageCounts       TextFile
sampleCumulativeCoverageProportions  TextFile
sampleIntervalStatistics             TextFile
sampleIntervalSummary                TextFile
sampleStatistics                     TextFile
sampleSummary                        TextFile
===================================  ========  ===============


Additional configuration (inputs)
---------------------------------

===============================  ========================  =================================  ==========  =====================================================================================================================
name                             type                      prefix                               position  documentation
===============================  ========================  =================================  ==========  =====================================================================================================================
bam                              IndexedBam                -I                                         10  Input file containing sequence  data (BAM or CRAM)
reference                        FastaWithIndexes          -R                                             Reference sequence file
outputPrefix                     String                    -o                                             An output file created by the walker. Will overwrite contents if file exists
intervals                        Optional<File>            -L                                             One or more genomic intervals over which to operate
excludeIntervals                 Optional<File>            --excludeIntervals                             One or more genomic intervals to exclude from processing
argFile                          Optional<File>            --arg_file                                     Reads arguments from the specified file
showFullBamList                  Optional<Boolean>         --showFullBamList                              Emit list of input BAM/CRAM files to log
read_buffer_size                 Optional<Integer>         --read_buffer_size                             Number of reads per SAM file to buffer in memory
read_filter                      Optional<Boolean>         --read_filter                                  Filters to apply to reads before analysis
disable_read_filter              Optional<Boolean>         --disable_read_filter                          Read filters to disable
interval_set_rule                Optional<String>          --interval_set_rule                            Set merging approach to use for combining interval inputs (UNION|INTERSECTION)
interval_merging                 Optional<String>          --interval_merging                             Set merging approach to use for combining interval inputs (UNION|INTERSECTION)
interval_padding                 Optional<Integer>         --interval_padding                             Amount of padding (in bp) to add to each interval
nonDeterministicRandomSeed       Optional<Boolean>         --nonDeterministicRandomSeed                   Use a non-deterministic random seed
maxRuntime                       Optional<String>          --maxRuntime                                   Unit of time used by maxRuntime (NANOSECONDS|MICROSECONDS|SECONDS|MINUTES|HOURS|DAYS)
downsampling_type                Optional<String>          --downsampling_type                            Type of read downsampling to employ at a given locus (NONE|ALL_READS|BY.sample)
downsample_to_fraction           Optional<Float>           --downsample_to_fraction                       Fraction of reads to downsample to Target coverage threshold for downsampling to coverage
baq                              Optional<String>          --baq                                          Type of BAQ calculation to apply in the engine (OFF|CALCULATE_AS_NECESSARY|RECALCULATE)
refactor_NDN_cigar_string        Optional<Boolean>         --refactor_NDN_cigar_string                    Reduce NDN elements in CIGAR string
fixMisencodedQuals               Optional<Boolean>         --fixMisencodedQuals                           Fix mis-encoded base quality scores
allowPotentiallyMisencodedQuals  Optional<Boolean>         --allowPotentiallyMisencodedQuals              Ignore warnings about base quality score encoding
useOriginalQualities             Optional<Boolean>         --useOriginalQualities                         Use the base quality scores from the OQ tag
defaultBaseQualities             Optional<Integer>         --defaultBaseQualities                         Assign a default base quality
performanceLog                   Optional<Filename>        --performanceLog                               Write GATK runtime performance log to this file
BQSR                             Optional<File>            --BQSR                                         Input covariates table file for on-the-fly base quality score recalibration
disable_indel_quals              Optional<Boolean>         --disable_indel_quals                          Disable printing of base insertion and deletion tags (with -BQSR)
emit_original_quals              Optional<Boolean>         --emit_original_quals                          Emit the OQ tag with the original base qualities (with -BQSR)
preserve_qscores_less_than       Optional<Integer>         --preserve_qscores_less_than                   Don't recalibrate bases with quality scores less than this threshold (with -BQSR)
countType                        Optional<String>          --countType                                    overlapping reads from the same  fragment be handled? (COUNT_READS|COUNT_FRAGMENTS|COUNT_FRAGMENTS_REQUIRE_SAME_BASE)
summaryCoverageThreshold         Optional<Array<Integer>>  -ct                                            Coverage threshold (in percent) for summarizing statistics
===============================  ========================  =================================  ==========  =====================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk3DepthOfCoverage {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File bam
       File bam_bai
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       String outputPrefix
       File? intervals
       File? excludeIntervals
       File? argFile
       Boolean? showFullBamList
       Int? read_buffer_size
       Boolean? read_filter
       Boolean? disable_read_filter
       String? interval_set_rule
       String? interval_merging
       Int? interval_padding
       Boolean? nonDeterministicRandomSeed
       String? maxRuntime
       String? downsampling_type
       Float? downsample_to_fraction
       String? baq
       Boolean? refactor_NDN_cigar_string
       Boolean? fixMisencodedQuals
       Boolean? allowPotentiallyMisencodedQuals
       Boolean? useOriginalQualities
       Int? defaultBaseQualities
       String? performanceLog
       File? BQSR
       Boolean? disable_indel_quals
       Boolean? emit_original_quals
       Int? preserve_qscores_less_than
       String? countType
       Array[Int]? summaryCoverageThreshold
     }
     command <<<
       set -e
       cp -f '~{bam_bai}' $(echo '~{bam}' | sed 's/\.[^.]*$//').bai
       java \
         -Xmx~{((select_first([runtime_memory, 4]) * 3) / 4)}G \
         -jar /usr/GenomeAnalysisTK.jar \
         -T DepthOfCoverage \
         -R '~{reference}' \
         -o '~{outputPrefix}' \
         ~{if defined(intervals) then ("-L '" + intervals + "'") else ""} \
         ~{if defined(excludeIntervals) then ("--excludeIntervals '" + excludeIntervals + "'") else ""} \
         ~{if defined(argFile) then ("--arg_file '" + argFile + "'") else ""} \
         ~{if (defined(showFullBamList) && select_first([showFullBamList])) then "--showFullBamList" else ""} \
         ~{if defined(read_buffer_size) then ("--read_buffer_size " + read_buffer_size) else ''} \
         ~{if (defined(read_filter) && select_first([read_filter])) then "--read_filter" else ""} \
         ~{if (defined(disable_read_filter) && select_first([disable_read_filter])) then "--disable_read_filter" else ""} \
         ~{if defined(interval_set_rule) then ("--interval_set_rule '" + interval_set_rule + "'") else ""} \
         ~{if defined(interval_merging) then ("--interval_merging '" + interval_merging + "'") else ""} \
         ~{if defined(interval_padding) then ("--interval_padding " + interval_padding) else ''} \
         ~{if (defined(nonDeterministicRandomSeed) && select_first([nonDeterministicRandomSeed])) then "--nonDeterministicRandomSeed" else ""} \
         ~{if defined(maxRuntime) then ("--maxRuntime '" + maxRuntime + "'") else ""} \
         ~{if defined(downsampling_type) then ("--downsampling_type '" + downsampling_type + "'") else ""} \
         ~{if defined(downsample_to_fraction) then ("--downsample_to_fraction " + downsample_to_fraction) else ''} \
         ~{if defined(baq) then ("--baq '" + baq + "'") else ""} \
         ~{if (defined(refactor_NDN_cigar_string) && select_first([refactor_NDN_cigar_string])) then "--refactor_NDN_cigar_string" else ""} \
         ~{if (defined(fixMisencodedQuals) && select_first([fixMisencodedQuals])) then "--fixMisencodedQuals" else ""} \
         ~{if (defined(allowPotentiallyMisencodedQuals) && select_first([allowPotentiallyMisencodedQuals])) then "--allowPotentiallyMisencodedQuals" else ""} \
         ~{if (defined(useOriginalQualities) && select_first([useOriginalQualities])) then "--useOriginalQualities" else ""} \
         ~{if defined(defaultBaseQualities) then ("--defaultBaseQualities " + defaultBaseQualities) else ''} \
         --performanceLog '~{select_first([performanceLog, "generated"])}' \
         ~{if defined(BQSR) then ("--BQSR '" + BQSR + "'") else ""} \
         ~{if (defined(disable_indel_quals) && select_first([disable_indel_quals])) then "--disable_indel_quals" else ""} \
         ~{if (defined(emit_original_quals) && select_first([emit_original_quals])) then "--emit_original_quals" else ""} \
         ~{if defined(preserve_qscores_less_than) then ("--preserve_qscores_less_than " + preserve_qscores_less_than) else ''} \
         ~{if defined(countType) then ("--countType '" + countType + "'") else ""} \
         ~{if (defined(summaryCoverageThreshold) && length(select_first([summaryCoverageThreshold])) > 0) then sep(" ", prefix("-ct ", select_first([summaryCoverageThreshold]))) else ""} \
         -I '~{bam}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk3:3.8-0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File sample = outputPrefix
       File sampleCumulativeCoverageCounts = (outputPrefix + ".sample_cumulative_coverage_counts")
       File sampleCumulativeCoverageProportions = (outputPrefix + ".sample_cumulative_coverage_proportions")
       File sampleIntervalStatistics = (outputPrefix + ".sample_interval_statistics")
       File sampleIntervalSummary = (outputPrefix + ".sample_interval_summary")
       File sampleStatistics = (outputPrefix + ".sample_statistics")
       File sampleSummary = (outputPrefix + ".sample_summary")
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: |-
     GATK3 DepthOfCoverage: Determine coverage at different levels of partitioning and aggregation.
   doc: |-
     Overview
     This tool processes a set of bam files to determine coverage at different levels of partitioning and aggregation. Coverage can be analyzed per locus, per interval, per gene, or in total; can be partitioned by sample, by read group, by technology, by center, or by library; and can be summarized by mean, median, quartiles, and/or percentage of bases covered to or beyond a threshold. Additionally, reads and bases can be filtered by mapping or base quality score.

     Input
     One or more bam files (with proper headers) to be analyzed for coverage statistics
     (Optional) A REFSEQ file to aggregate coverage to the gene level (for information about creating the REFSEQ Rod, please consult the online documentation)
     Output
     Tables pertaining to different coverage summaries. Suffix on the table files declares the contents:

     no suffix: per locus coverage
     _summary: total, mean, median, quartiles, and threshold proportions, aggregated over all bases
     _statistics: coverage histograms (# locus with X coverage), aggregated over all bases
     _interval_summary: total, mean, median, quartiles, and threshold proportions, aggregated per interval
     _interval_statistics: 2x2 table of # of intervals covered to >= X depth in >=Y samples
     _gene_summary: total, mean, median, quartiles, and threshold proportions, aggregated per gene
     _gene_statistics: 2x2 table of # of genes covered to >= X depth in >= Y samples
     _cumulative_coverage_counts: coverage histograms (# locus with >= X coverage), aggregated over all bases
     _cumulative_coverage_proportions: proprotions of loci with >= X coverage, aggregated over all bases

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk3:3.8-0

   inputs:
   - id: bam
     label: bam
     doc: Input file containing sequence  data (BAM or CRAM)
     type: File
     secondaryFiles:
     - |-
       ${

               function resolveSecondary(base, secPattern) {
                 if (secPattern[0] == "^") {
                   var spl = base.split(".");
                   var endIndex = spl.length > 1 ? spl.length - 1 : 1;
                   return resolveSecondary(spl.slice(undefined, endIndex).join("."), secPattern.slice(1));
                 }
                 return base + secPattern
               }

               return [
                       {
                           location: resolveSecondary(self.location, "^.bai"),
                           basename: resolveSecondary(self.basename, ".bai"),
                           class: "File",
                       }
               ];

       }
     inputBinding:
       prefix: -I
       position: 10
   - id: reference
     label: reference
     doc: Reference sequence file
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
       prefix: -R
   - id: outputPrefix
     label: outputPrefix
     doc: An output file created by the walker. Will overwrite contents if file exists
     type: string
     inputBinding:
       prefix: -o
   - id: intervals
     label: intervals
     doc: One or more genomic intervals over which to operate
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -L
   - id: excludeIntervals
     label: excludeIntervals
     doc: One or more genomic intervals to exclude from processing
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --excludeIntervals
   - id: argFile
     label: argFile
     doc: Reads arguments from the specified file
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --arg_file
   - id: showFullBamList
     label: showFullBamList
     doc: Emit list of input BAM/CRAM files to log
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --showFullBamList
   - id: read_buffer_size
     label: read_buffer_size
     doc: Number of reads per SAM file to buffer in memory
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --read_buffer_size
   - id: read_filter
     label: read_filter
     doc: Filters to apply to reads before analysis
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --read_filter
   - id: disable_read_filter
     label: disable_read_filter
     doc: Read filters to disable
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --disable_read_filter
   - id: interval_set_rule
     label: interval_set_rule
     doc: Set merging approach to use for combining interval inputs (UNION|INTERSECTION)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --interval_set_rule
   - id: interval_merging
     label: interval_merging
     doc: Set merging approach to use for combining interval inputs (UNION|INTERSECTION)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --interval_merging
   - id: interval_padding
     label: interval_padding
     doc: Amount of padding (in bp) to add to each interval
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --interval_padding
   - id: nonDeterministicRandomSeed
     label: nonDeterministicRandomSeed
     doc: Use a non-deterministic random seed
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --nonDeterministicRandomSeed
   - id: maxRuntime
     label: maxRuntime
     doc: |-
       Unit of time used by maxRuntime (NANOSECONDS|MICROSECONDS|SECONDS|MINUTES|HOURS|DAYS)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --maxRuntime
   - id: downsampling_type
     label: downsampling_type
     doc: Type of read downsampling to employ at a given locus (NONE|ALL_READS|BY.sample)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --downsampling_type
   - id: downsample_to_fraction
     label: downsample_to_fraction
     doc: |-
       Fraction of reads to downsample to Target coverage threshold for downsampling to coverage
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --downsample_to_fraction
   - id: baq
     label: baq
     doc: |-
       Type of BAQ calculation to apply in the engine (OFF|CALCULATE_AS_NECESSARY|RECALCULATE)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --baq
   - id: refactor_NDN_cigar_string
     label: refactor_NDN_cigar_string
     doc: Reduce NDN elements in CIGAR string
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --refactor_NDN_cigar_string
   - id: fixMisencodedQuals
     label: fixMisencodedQuals
     doc: Fix mis-encoded base quality scores
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --fixMisencodedQuals
   - id: allowPotentiallyMisencodedQuals
     label: allowPotentiallyMisencodedQuals
     doc: Ignore warnings about base quality score encoding
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --allowPotentiallyMisencodedQuals
   - id: useOriginalQualities
     label: useOriginalQualities
     doc: Use the base quality scores from the OQ tag
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --useOriginalQualities
   - id: defaultBaseQualities
     label: defaultBaseQualities
     doc: Assign a default base quality
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --defaultBaseQualities
   - id: performanceLog
     label: performanceLog
     doc: Write GATK runtime performance log to this file
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --performanceLog
   - id: BQSR
     label: BQSR
     doc: Input covariates table file for on-the-fly base quality score recalibration
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --BQSR
   - id: disable_indel_quals
     label: disable_indel_quals
     doc: Disable printing of base insertion and deletion tags (with -BQSR)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --disable_indel_quals
   - id: emit_original_quals
     label: emit_original_quals
     doc: Emit the OQ tag with the original base qualities (with -BQSR)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --emit_original_quals
   - id: preserve_qscores_less_than
     label: preserve_qscores_less_than
     doc: |-
       Don't recalibrate bases with quality scores less than this threshold (with -BQSR)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --preserve_qscores_less_than
   - id: countType
     label: countType
     doc: |-
       overlapping reads from the same  fragment be handled? (COUNT_READS|COUNT_FRAGMENTS|COUNT_FRAGMENTS_REQUIRE_SAME_BASE)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --countType
   - id: summaryCoverageThreshold
     label: summaryCoverageThreshold
     doc: Coverage threshold (in percent) for summarizing statistics
     type:
     - type: array
       inputBinding:
         prefix: -ct
       items: int
     - 'null'
     inputBinding: {}

   outputs:
   - id: sample
     label: sample
     doc: ''
     type: File
     outputBinding:
       glob: $(inputs.outputPrefix)
       loadContents: false
   - id: sampleCumulativeCoverageCounts
     label: sampleCumulativeCoverageCounts
     doc: ''
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".sample_cumulative_coverage_counts"))
       outputEval: $((inputs.outputPrefix.basename + ".sample_cumulative_coverage_counts"))
       loadContents: false
   - id: sampleCumulativeCoverageProportions
     label: sampleCumulativeCoverageProportions
     doc: ''
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".sample_cumulative_coverage_proportions"))
       outputEval: $((inputs.outputPrefix.basename + ".sample_cumulative_coverage_proportions"))
       loadContents: false
   - id: sampleIntervalStatistics
     label: sampleIntervalStatistics
     doc: ''
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".sample_interval_statistics"))
       outputEval: $((inputs.outputPrefix.basename + ".sample_interval_statistics"))
       loadContents: false
   - id: sampleIntervalSummary
     label: sampleIntervalSummary
     doc: ''
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".sample_interval_summary"))
       outputEval: $((inputs.outputPrefix.basename + ".sample_interval_summary"))
       loadContents: false
   - id: sampleStatistics
     label: sampleStatistics
     doc: ''
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".sample_statistics"))
       outputEval: $((inputs.outputPrefix.basename + ".sample_statistics"))
       loadContents: false
   - id: sampleSummary
     label: sampleSummary
     doc: ''
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".sample_summary"))
       outputEval: $((inputs.outputPrefix.basename + ".sample_summary"))
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - java
   arguments:
   - position: -3
     valueFrom: |-
       $("-Xmx{memory}G".replace(/\{memory\}/g, (([inputs.runtime_memory, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)))
     shellQuote: false
   - position: -2
     valueFrom: -jar /usr/GenomeAnalysisTK.jar
     shellQuote: false
   - position: -1
     valueFrom: -T DepthOfCoverage
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk3DepthOfCoverage


