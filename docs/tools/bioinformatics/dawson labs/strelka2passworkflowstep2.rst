:orphan:

Strelka 2Pass analysis step 2
=========================================================

*1 contributor · 1 version*

This is the second step for joint somatic variant calling
        based on a 2pass analysis common in RNASeq.

        It runs strelka2 again with the variants found in all of the other samples as input to be forced to genotype these.

        It also normalises and indexes the output vcfs


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.dawson.workflows.strelka2passanalysisstep2 import Strelka2PassWorkflowStep2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "strelka2passworkflowstep2_step",
           Strelka2PassWorkflowStep2(
               normalBam=None,
               tumorBam=None,
               reference=None,
               indelCandidates=None,
               strelkaSNVs=None,
           )
       )
       wf.output("indels", source=strelka2passworkflowstep2_step.indels)
       wf.output("snvs", source=strelka2passworkflowstep2_step.snvs)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Strelka2PassWorkflowStep2:

.. code-block:: bash

   # user inputs
   janis inputs Strelka2PassWorkflowStep2 > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       indelCandidates:
       - indelCandidates_0.vcf.gz
       - indelCandidates_1.vcf.gz
       normalBam: normalBam.cram
       reference: reference.fasta
       strelkaSNVs:
       - strelkaSNVs_0.vcf.gz
       - strelkaSNVs_1.vcf.gz
       tumorBam: tumorBam.cram




5. Run Strelka2PassWorkflowStep2 with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Strelka2PassWorkflowStep2





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``Strelka2PassWorkflowStep2``
:URL: *No URL to the documentation was provided*
:Versions: 0.1
:Authors: Sebastian Hollizeck
:Citations: 
:Created: 2019-10-11
:Updated: 2019-10-11



Outputs
-----------

======  ====================  ===============
name    type                  documentation
======  ====================  ===============
indels  CompressedIndexedVCF
snvs    CompressedIndexedVCF
======  ====================  ===============


Embedded Tools
***************

===================  ===============================
Strelka (Somatic)    ``strelka_somatic_cram/2.9.10``
BCFTools: Normalize  ``bcftoolsNorm/v1.9``
BCFTools: Index      ``bcftoolsIndex/v1.9``
===================  ===============================



Additional configuration (inputs)
---------------------------------

===============  ===========================  ===============
name             type                         documentation
===============  ===========================  ===============
normalBam        CramPair
tumorBam         CramPair
reference        FastaWithIndexes
indelCandidates  Array<CompressedIndexedVCF>
strelkaSNVs      Array<CompressedIndexedVCF>
callRegions      Optional<BedTABIX>
exome            Optional<Boolean>
===============  ===========================  ===============


