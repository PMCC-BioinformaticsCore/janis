:orphan:

GATK3 DepthOfCoverage: Determine coverage at different levels of partitioning and aggregation.
=====================================================================================================================

*1 contributor · 1 version*

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

       from janis_bioinformatics.tools.gatk3.depthofcoverage.versions import GATK3DepthOfCoverage_3_8_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk3depthofcoverage_step",
           GATK3DepthOfCoverage_3_8_1(
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
:Versions: 3.8-1
:Container: broadinstitute/gatk3:3.8-1
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
