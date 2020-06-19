:orphan:

Strelka 2Pass analysis step1
========================================================

*1 contributor · 1 version*

This is the first step for joint somatic variant calling
        based on a 2pass analysis common in RNASeq.

        It runs manta and strelka on the bams as is best practice
        for somatic variant calling with strelka2

        It also normalises and indexes the output vcfs


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.dawson.workflows.strelka2passanalysisstep1 import Strelka2PassWorkflowStep1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "strelka2passworkflowstep1_step",
           Strelka2PassWorkflowStep1(
               normalBam=None,
               tumorBam=None,
               reference=None,
           )
       )
       wf.output("diploid", source=strelka2passworkflowstep1_step.diploid)
       wf.output("candIndels", source=strelka2passworkflowstep1_step.candIndels)
       wf.output("indels", source=strelka2passworkflowstep1_step.indels)
       wf.output("snvs", source=strelka2passworkflowstep1_step.snvs)
       wf.output("somaticSVs", source=strelka2passworkflowstep1_step.somaticSVs)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Strelka2PassWorkflowStep1:

.. code-block:: bash

   # user inputs
   janis inputs Strelka2PassWorkflowStep1 > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       normalBam: normalBam.cram
       reference: reference.fasta
       tumorBam: tumorBam.cram




5. Run Strelka2PassWorkflowStep1 with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Strelka2PassWorkflowStep1





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``Strelka2PassWorkflowStep1``
:URL: *No URL to the documentation was provided*
:Versions: 0.1
:Authors: Sebastian Hollizeck
:Citations: 
:Created: 2019-10-11
:Updated: 2019-10-11



Outputs
-----------

==========  ====================  ===============
name        type                  documentation
==========  ====================  ===============
diploid     CompressedIndexedVCF
candIndels  CompressedIndexedVCF
indels      CompressedIndexedVCF
snvs        CompressedIndexedVCF
somaticSVs  CompressedIndexedVCF
==========  ====================  ===============


Embedded Tools
***************

===================  ===============================
Manta                ``manta_cram/1.5.0``
Strelka (Somatic)    ``strelka_somatic_cram/2.9.10``
BCFTools: Normalize  ``bcftoolsNorm/v1.9``
BCFTools: Index      ``bcftoolsIndex/v1.9``
===================  ===============================



Additional configuration (inputs)
---------------------------------

===========  ==================  ===============
name         type                documentation
===========  ==================  ===============
normalBam    CramPair
tumorBam     CramPair
reference    FastaWithIndexes
callRegions  Optional<BedTABIX>
exome        Optional<Boolean>
===========  ==================  ===============


