:orphan:

GATK4: Generate coverage summary information for reads data
==================================================================================

*1 contributor · 1 version*

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
