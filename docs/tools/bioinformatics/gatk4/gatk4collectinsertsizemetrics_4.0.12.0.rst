:orphan:

GATK4: CollectInsertSizeMetrics
===============================================================

*1 contributor · 4 versions*

Provides useful metrics for validating library construction including the insert size distribution and read orientation of paired-end libraries


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.collectinsertsizemetrics.versions import Gatk4CollectInsertSizeMetrics_4_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4collectinsertsizemetrics_step",
           Gatk4CollectInsertSizeMetrics_4_0(
               bam=None,
           )
       )
       wf.output("out", source=gatk4collectinsertsizemetrics_step.out)
       wf.output("outHistogram", source=gatk4collectinsertsizemetrics_step.outHistogram)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4CollectInsertSizeMetrics:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4CollectInsertSizeMetrics > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam




5. Run Gatk4CollectInsertSizeMetrics with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4CollectInsertSizeMetrics





Information
------------


:ID: ``Gatk4CollectInsertSizeMetrics``
:URL: `https://gatk.broadinstitute.org/hc/en-us/articles/360036715591-CollectInsertSizeMetrics-Picard- <https://gatk.broadinstitute.org/hc/en-us/articles/360036715591-CollectInsertSizeMetrics-Picard->`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.0.12.0
:Authors: Jiaan Yu
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2020-02-17
:Updated: 2020-02-17



Outputs
-----------

============  ========  ===============
name          type      documentation
============  ========  ===============
out           TextFile
outHistogram  File
============  ========  ===============



Additional configuration (inputs)
---------------------------------

=======================  =======================  ===========================  ==========  ============================================================================================================================================================================================================================================================================================================================
name                     type                     prefix                         position  documentation
=======================  =======================  ===========================  ==========  ============================================================================================================================================================================================================================================================================================================================
bam                      IndexedBam               -I                                   10  Input SAM or BAM file.  Required.
javaOptions              Optional<Array<String>>
compression_level        Optional<Integer>                                                 Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
outputFilename           Optional<Filename>       -O                                       File to write the output to.  Required.
outputHistogram          Optional<Filename>       -H                                       File to write insert size Histogram chart to.  Required.
argumentsFile            Optional<Array<File>>    --arguments_file                     10  read one or more arguments files and add them to the command line
assumeSorted             Optional<Boolean>        --ASSUME_SORTED                      11  If true (default), then the sort order in the header file will be ignored.  Default value: true. Possible values: {true, false}
deviations               Optional<Double>         --DEVIATIONS                         11  Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. This is done because insert size data typically includes enough anomalous values from chimeras and other artifacts to make the mean and sd grossly misleading regarding the real distribution.  Default value: 10.0.
histogramWidth           Optional<Integer>        --HISTOGRAM_WIDTH                    11  Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. Also, when calculating mean and standard deviation, only bins <= Histogram_WIDTH will be included.  Default value: null.
includeDuplicates        Optional<Boolean>        --INCLUDE_DUPLICATES                 11  If true, also include reads marked as duplicates in the insert size histogram.  Default value: false. Possible values: {true, false}
metricAccumulationLevel  Optional<String>         --METRIC_ACCUMULATION_LEVEL          11  The level(s) at  which to accumulate metrics.    This argument may be specified 0 or more times. Default value: [ALL_READS]. Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} .
minimumPCT               Optional<Float>          --MINIMUM_PCT                        11  When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this percentage of overall reads. (Range: 0 to 1).  Default value: 0.05.
stopAfter                Optional<Integer>        --STOP_AFTER                         11  Stop after  processing N reads, mainly for debugging.  Default value: 0.
version                  Optional<Boolean>        --version                            11  display the version number for this tool Default value: false. Possible values: {true, false}
showHidden               Optional<Boolean>        --showHidden                         11  display hidden  arguments  Default  value: false.  Possible values: {true, false}
=======================  =======================  ===========================  ==========  ============================================================================================================================================================================================================================================================================================================================
