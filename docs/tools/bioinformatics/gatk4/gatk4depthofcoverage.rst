:orphan:

GATK4: Generate coverage summary information for reads data
==================================================================================

``Gatk4DepthOfCoverage`` · *1 contributor · 1 version*

Generate coverage summary information for reads data

Category Coverage Analysis
Overview
Assess sequence coverage by a wide array of metrics, partitioned by sample, read group, or library
This tool processes a set of bam files to determine coverage at different levels of partitioning and aggregation. Coverage can be analyzed per locus, per interval, per gene, or in total; can be partitioned by sample, by read group, by technology, by center, or by library; and can be summarized by mean, median, quartiles, and/or percentage of bases covered to or beyond a threshold. Additionally, reads and bases can be filtered by mapping or base quality score.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.depthofcoverage.versions import Gatk4DepthOfCoverage_4_1_6

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4depthofcoverage_step",
           Gatk4DepthOfCoverage_4_1_6(
               bam=None,
               reference=None,
               outputPrefix=None,
               intervals=None,
           )
       )
       wf.output("out_sample", source=gatk4depthofcoverage_step.out_sample)
       wf.output("out_sampleCumulativeCoverageCounts", source=gatk4depthofcoverage_step.out_sampleCumulativeCoverageCounts)
       wf.output("out_sampleCumulativeCoverageProportions", source=gatk4depthofcoverage_step.out_sampleCumulativeCoverageProportions)
       wf.output("out_sampleIntervalStatistics", source=gatk4depthofcoverage_step.out_sampleIntervalStatistics)
       wf.output("out_sampleIntervalSummary", source=gatk4depthofcoverage_step.out_sampleIntervalSummary)
       wf.output("out_sampleStatistics", source=gatk4depthofcoverage_step.out_sampleStatistics)
       wf.output("out_sampleSummary", source=gatk4depthofcoverage_step.out_sampleSummary)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4DepthOfCoverage:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4DepthOfCoverage > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       intervals:
       - intervals_0.bed
       - intervals_1.bed
       outputPrefix: <value>
       reference: reference.fasta




5. Run Gatk4DepthOfCoverage with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4DepthOfCoverage





Information
------------

:ID: ``Gatk4DepthOfCoverage``
:URL: `https://gatk.broadinstitute.org/hc/en-us/articles/360041851491-DepthOfCoverage-BETA- <https://gatk.broadinstitute.org/hc/en-us/articles/360041851491-DepthOfCoverage-BETA->`_
:Versions: 4.1.6.0
:Container: broadinstitute/gatk:4.1.6.0
:Authors: Jiaan Yu
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2020-07-10
:Updated: 2020-07-10


Outputs
-----------

=======================================  ==================  ====================================================================================
name                                     type                documentation
=======================================  ==================  ====================================================================================
out_sample                               Optional<TextFile>  per locus coverage
out_sampleCumulativeCoverageCounts       TextFile            coverage histograms (# locus with >= X coverage), aggregated over all bases
out_sampleCumulativeCoverageProportions  TextFile            proprotions of loci with >= X coverage, aggregated over all bases
out_sampleIntervalStatistics             TextFile            total, mean, median, quartiles, and threshold proportions, aggregated per interval
out_sampleIntervalSummary                TextFile            2x2 table of # of intervals covered to >= X depth in >=Y samples
out_sampleStatistics                     TextFile            coverage histograms (# locus with X coverage), aggregated over all bases
out_sampleSummary                        TextFile            total, mean, median, quartiles, and threshold proportions, aggregated over all bases
=======================================  ==================  ====================================================================================


Additional configuration (inputs)
---------------------------------

======================================  ========================  ==============================================  ==========  =====================================================================================================================
name                                    type                      prefix                                          position    documentation
======================================  ========================  ==============================================  ==========  =====================================================================================================================
bam                                     IndexedBam                -I                                                          The SAM/BAM/CRAM file containing reads.
reference                               FastaWithIndexes          -R                                                          Reference sequence
outputPrefix                            String                    -O                                                          An output file created by the walker. Will overwrite contents if file exists
intervals                               Array<bed>                --intervals                                                 -L (BASE) One or more genomic intervals over which to operate
javaOptions                             Optional<Array<String>>
compression_level                       Optional<Integer>                                                                     Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
countType                               Optional<String>          --count-type                                                overlapping reads from the same  fragment be handled? (COUNT_READS|COUNT_FRAGMENTS|COUNT_FRAGMENTS_REQUIRE_SAME_BASE)
summaryCoverageThreshold                Optional<Array<Integer>>  --summary-coverage-threshold                                Coverage threshold (in percent) for summarizing statistics
omitDepthOutputAtEachBase               Optional<Boolean>         --omit-depth-output-at-each-base                            Do not output depth of coverage at each base
omitGenesNotEntirelyCoveredByTraversal  Optional<Boolean>         --omit-genes-not-entirely-covered-by-traversal              Do not output gene summary if it was not completely covered by traversal intervals
omitIntervalStatistics                  Optional<Boolean>         --omit-interval-statistics                                  Do not calculate per-interval statistics
omitLocusTable                          Optional<Boolean>         --omit-locus-table                                          Do not calculate per-sample per-depth counts of loci
omitPerSampleStatistics                 Optional<Boolean>         --omit-per-sample-statistics                                Do not output the summary files per-sample
======================================  ========================  ==============================================  ==========  =====================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4DepthOfCoverage {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
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
       Array[File] intervals
       String? countType
       Array[Int]? summaryCoverageThreshold
       Boolean? omitDepthOutputAtEachBase
       Boolean? omitGenesNotEntirelyCoveredByTraversal
       Boolean? omitIntervalStatistics
       Boolean? omitLocusTable
       Boolean? omitPerSampleStatistics
     }
     command <<<
       set -e
       cp -f '~{bam_bai}' $(echo '~{bam}' | sed 's/\.[^.]*$//').bai
       gatk DepthOfCoverage \
         --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         -I '~{bam}' \
         -R '~{reference}' \
         -O '~{outputPrefix}' \
         ~{if length(intervals) > 0 then "--intervals '" + sep("' --intervals '", intervals) + "'" else ""} \
         ~{if defined(countType) then ("--count-type '" + countType + "'") else ""} \
         ~{if (defined(summaryCoverageThreshold) && length(select_first([summaryCoverageThreshold])) > 0) then sep(" ", prefix("--summary-coverage-threshold ", select_first([summaryCoverageThreshold]))) else ""} \
         ~{if (defined(omitDepthOutputAtEachBase) && select_first([omitDepthOutputAtEachBase])) then "--omit-depth-output-at-each-base" else ""} \
         ~{if (defined(omitGenesNotEntirelyCoveredByTraversal) && select_first([omitGenesNotEntirelyCoveredByTraversal])) then "--omit-genes-not-entirely-covered-by-traversal" else ""} \
         ~{if (defined(omitIntervalStatistics) && select_first([omitIntervalStatistics])) then "--omit-interval-statistics" else ""} \
         ~{if (defined(omitLocusTable) && select_first([omitLocusTable])) then "--omit-locus-table" else ""} \
         ~{if (defined(omitPerSampleStatistics) && select_first([omitPerSampleStatistics])) then "--omit-per-sample-statistics" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.6.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File? out_sample = outputPrefix
       File out_sampleCumulativeCoverageCounts = (outputPrefix + ".sample_cumulative_coverage_counts")
       File out_sampleCumulativeCoverageProportions = (outputPrefix + ".sample_cumulative_coverage_proportions")
       File out_sampleIntervalStatistics = (outputPrefix + ".sample_interval_statistics")
       File out_sampleIntervalSummary = (outputPrefix + ".sample_interval_summary")
       File out_sampleStatistics = (outputPrefix + ".sample_statistics")
       File out_sampleSummary = (outputPrefix + ".sample_summary")
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: Generate coverage summary information for reads data'
   doc: |-
     Generate coverage summary information for reads data

     Category Coverage Analysis
     Overview
     Assess sequence coverage by a wide array of metrics, partitioned by sample, read group, or library
     This tool processes a set of bam files to determine coverage at different levels of partitioning and aggregation. Coverage can be analyzed per locus, per interval, per gene, or in total; can be partitioned by sample, by read group, by technology, by center, or by library; and can be summarized by mean, median, quartiles, and/or percentage of bases covered to or beyond a threshold. Additionally, reads and bases can be filtered by mapping or base quality score.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.6.0

   inputs:
   - id: javaOptions
     label: javaOptions
     type:
     - type: array
       items: string
     - 'null'
   - id: compression_level
     label: compression_level
     doc: |-
       Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
     type:
     - int
     - 'null'
   - id: bam
     label: bam
     doc: The SAM/BAM/CRAM file containing reads.
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
   - id: reference
     label: reference
     doc: Reference sequence
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
       prefix: -O
   - id: intervals
     label: intervals
     doc: -L (BASE) One or more genomic intervals over which to operate
     type:
       type: array
       inputBinding:
         prefix: --intervals
       items: File
     inputBinding: {}
   - id: countType
     label: countType
     doc: |-
       overlapping reads from the same  fragment be handled? (COUNT_READS|COUNT_FRAGMENTS|COUNT_FRAGMENTS_REQUIRE_SAME_BASE)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --count-type
   - id: summaryCoverageThreshold
     label: summaryCoverageThreshold
     doc: Coverage threshold (in percent) for summarizing statistics
     type:
     - type: array
       inputBinding:
         prefix: --summary-coverage-threshold
       items: int
     - 'null'
     inputBinding: {}
   - id: omitDepthOutputAtEachBase
     label: omitDepthOutputAtEachBase
     doc: Do not output depth of coverage at each base
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --omit-depth-output-at-each-base
   - id: omitGenesNotEntirelyCoveredByTraversal
     label: omitGenesNotEntirelyCoveredByTraversal
     doc: |-
       Do not output gene summary if it was not completely covered by traversal intervals
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --omit-genes-not-entirely-covered-by-traversal
   - id: omitIntervalStatistics
     label: omitIntervalStatistics
     doc: Do not calculate per-interval statistics
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --omit-interval-statistics
   - id: omitLocusTable
     label: omitLocusTable
     doc: Do not calculate per-sample per-depth counts of loci
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --omit-locus-table
   - id: omitPerSampleStatistics
     label: omitPerSampleStatistics
     doc: Do not output the summary files per-sample
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --omit-per-sample-statistics

   outputs:
   - id: out_sample
     label: out_sample
     doc: per locus coverage
     type:
     - File
     - 'null'
     outputBinding:
       glob: $(inputs.outputPrefix)
       loadContents: false
   - id: out_sampleCumulativeCoverageCounts
     label: out_sampleCumulativeCoverageCounts
     doc: coverage histograms (# locus with >= X coverage), aggregated over all bases
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".sample_cumulative_coverage_counts"))
       outputEval: $((inputs.outputPrefix.basename + ".sample_cumulative_coverage_counts"))
       loadContents: false
   - id: out_sampleCumulativeCoverageProportions
     label: out_sampleCumulativeCoverageProportions
     doc: proprotions of loci with >= X coverage, aggregated over all bases
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".sample_cumulative_coverage_proportions"))
       outputEval: $((inputs.outputPrefix.basename + ".sample_cumulative_coverage_proportions"))
       loadContents: false
   - id: out_sampleIntervalStatistics
     label: out_sampleIntervalStatistics
     doc: |-
       total, mean, median, quartiles, and threshold proportions, aggregated per interval
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".sample_interval_statistics"))
       outputEval: $((inputs.outputPrefix.basename + ".sample_interval_statistics"))
       loadContents: false
   - id: out_sampleIntervalSummary
     label: out_sampleIntervalSummary
     doc: '2x2 table of # of intervals covered to >= X depth in >=Y samples'
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".sample_interval_summary"))
       outputEval: $((inputs.outputPrefix.basename + ".sample_interval_summary"))
       loadContents: false
   - id: out_sampleStatistics
     label: out_sampleStatistics
     doc: coverage histograms (# locus with X coverage), aggregated over all bases
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".sample_statistics"))
       outputEval: $((inputs.outputPrefix.basename + ".sample_statistics"))
       loadContents: false
   - id: out_sampleSummary
     label: out_sampleSummary
     doc: |-
       total, mean, median, quartiles, and threshold proportions, aggregated over all bases
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".sample_summary"))
       outputEval: $((inputs.outputPrefix.basename + ".sample_summary"))
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - DepthOfCoverage
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4DepthOfCoverage


