:orphan:

Performance summary workflow (targeted bed)
========================================================================

*1 contributor · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.performanceSummaryTargetedWorkflow import PerformanceSummaryTargeted_0_1_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "performancesummarytargeted_step",
           PerformanceSummaryTargeted_0_1_0(
               bam=None,
               genecoverage_bed=None,
               region_bed=None,
               sample_name=None,
               genome_file=None,
           )
       )
       wf.output("out", source=performancesummarytargeted_step.out)
       wf.output("geneFileOut", source=performancesummarytargeted_step.geneFileOut)
       wf.output("regionFileOut", source=performancesummarytargeted_step.regionFileOut)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for PerformanceSummaryTargeted:

.. code-block:: bash

   # user inputs
   janis inputs PerformanceSummaryTargeted > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       genecoverage_bed: genecoverage_bed.bed
       genome_file: genome_file.txt
       region_bed: region_bed.bed
       sample_name: <value>




5. Run PerformanceSummaryTargeted with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       PerformanceSummaryTargeted





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``PerformanceSummaryTargeted``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Authors: Jiaan Yu
:Citations: 
:Created: None
:Updated: None



Outputs
-----------

=============  ========  ===============
name           type      documentation
=============  ========  ===============
out            csv
geneFileOut    TextFile
regionFileOut  TextFile
=============  ========  ===============


Embedded Tools
***************

===============================  =========================================
GATK4: CollectInsertSizeMetrics  ``Gatk4CollectInsertSizeMetrics/4.1.3.0``
SamTools: Flagstat               ``SamToolsFlagstat/1.9.0``
SamTools: View                   ``SamToolsView/1.9.0``
BEDTools: intersectBed           ``bedtoolsintersectBed/v2.29.2``
BEDTools: coverageBed            ``bedtoolsCoverageBed/v2.29.2``
Performance Summary              ``performanceSummary/0.0.7``
Gene Coverage Per Sample         ``geneCoveragePerSample/0.0.8``
===============================  =========================================



Additional configuration (inputs)
---------------------------------

=============================================  =================  ==========================================================================================================================================================================================================================
name                                           type               documentation
=============================================  =================  ==========================================================================================================================================================================================================================
bam                                            IndexedBam
genecoverage_bed                               bed
region_bed                                     bed
sample_name                                    String
genome_file                                    TextFile
samtoolsview_doNotOutputAlignmentsWithBitsSet  Optional<String>   Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].
bedtoolsintersectbed_sorted                    Optional<Boolean>  Use the 'chromsweep' algorithm for sorted (-k1,1 -k2,2n) input.
bedtoolscoveragebed_sorted                     Optional<Boolean>  Use the 'chromsweep' algorithm for sorted (-k1,1 -k2,2n) input.
bedtoolscoveragebed_histogram                  Optional<Boolean>  Report a histogram of coverage for each feature in A as well as a summary histogram for _all_ features in A. Output (tab delimited) after each feature in A: 1) depth 2) # bases at depth 3) size of A 4) % of A at depth.
bedtoolscoverage_sorted                        Optional<Boolean>  Use the 'chromsweep' algorithm for sorted (-k1,1 -k2,2n) input.
bedtoolscoverage_histogram                     Optional<Boolean>  Report a histogram of coverage for each feature in A as well as a summary histogram for _all_ features in A. Output (tab delimited) after each feature in A: 1) depth 2) # bases at depth 3) size of A 4) % of A at depth.
=============================================  =================  ==========================================================================================================================================================================================================================


