:orphan:

Annotate GATK3 DepthOfCoverage Workflow
=================================================================

*1 contributor · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.annotateDepthOfCoverageWorkflow import AnnotateDepthOfCoverage_0_1_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "annotatedepthofcoverage_step",
           AnnotateDepthOfCoverage_0_1_0(
               bam=None,
               bed=None,
               reference=None,
               sample_name=None,
           )
       )
       wf.output("out", source=annotatedepthofcoverage_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for AnnotateDepthOfCoverage:

.. code-block:: bash

   # user inputs
   janis inputs AnnotateDepthOfCoverage > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       bed: bed.bed
       reference: reference.fasta
       sample_name: <value>




5. Run AnnotateDepthOfCoverage with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       AnnotateDepthOfCoverage





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``AnnotateDepthOfCoverage``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Authors: Jiaan Yu
:Citations: 
:Created: None
:Updated: None



Outputs
-----------

======  ========  ===============
name    type      documentation
======  ========  ===============
out     TextFile
======  ========  ===============


Embedded Tools
***************

==============================================================================================  =================================
GATK3 DepthOfCoverage: Determine coverage at different levels of partitioning and aggregation.  ``Gatk3DepthOfCoverage/3.8-1``
Add Sym to DepthOfCoverage                                                                      ``addSymToDepthOfCoverage/0.0.7``
==============================================================================================  =================================



Additional configuration (inputs)
---------------------------------

=============================================  ========================  =====================================================================================================================
name                                           type                      documentation
=============================================  ========================  =====================================================================================================================
bam                                            IndexedBam
bed                                            bed
reference                                      FastaWithIndexes
sample_name                                    String
gatk3depthofcoverage_countType                 Optional<String>          overlapping reads from the same  fragment be handled? (COUNT_READS|COUNT_FRAGMENTS|COUNT_FRAGMENTS_REQUIRE_SAME_BASE)
gatk3depthofcoverage_summaryCoverageThreshold  Optional<Array<Integer>>  Coverage threshold (in percent) for summarizing statistics
=============================================  ========================  =====================================================================================================================


