:orphan:

Performance summary workflow (whole genome)
======================================================================

*1 contributor · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.performanceSummaryGenomeWorkflow import PerformanceSummaryGenome_0_1_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "performancesummarygenome_step",
           PerformanceSummaryGenome_0_1_0(
               bam=None,
               sample_name=None,
               genome_file=None,
           )
       )
       wf.output("performanceSummaryOut", source=performancesummarygenome_step.performanceSummaryOut)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for PerformanceSummaryGenome:

.. code-block:: bash

   # user inputs
   janis inputs PerformanceSummaryGenome > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       genome_file: genome_file.txt
       sample_name: <value>




5. Run PerformanceSummaryGenome with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       PerformanceSummaryGenome





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``PerformanceSummaryGenome``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Authors: Jiaan Yu
:Citations: 
:Created: None
:Updated: None



Outputs
-----------

=====================  ======  ===============
name                   type    documentation
=====================  ======  ===============
performanceSummaryOut  csv
=====================  ======  ===============


Embedded Tools
***************

===============================  =========================================
GATK4: CollectInsertSizeMetrics  ``Gatk4CollectInsertSizeMetrics/4.1.3.0``
SamTools: Flagstat               ``SamToolsFlagstat/1.9.0``
SamTools: View                   ``SamToolsView/1.9.0``
BEDTools: genomeCoverageBed      ``bedtoolsgenomeCoverageBed/v2.29.2``
Performance Summary              ``performanceSummary/0.0.7``
===============================  =========================================



Additional configuration (inputs)
---------------------------------

=============================================  =================  ==============================================================================================================================================================================================================
name                                           type               documentation
=============================================  =================  ==============================================================================================================================================================================================================
bam                                            IndexedBam
sample_name                                    String
genome_file                                    TextFile
samtoolsview_doNotOutputAlignmentsWithBitsSet  Optional<String>   Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].
performancesummary_genome                      Optional<Boolean>  calculate statistics for whole genome data.--target_flagstat must not be speicified
=============================================  =================  ==============================================================================================================================================================================================================


