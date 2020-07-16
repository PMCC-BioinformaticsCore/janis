:orphan:

Mutect2 joint somatic variant calling workflow
============================================================================

*0 contributors · 1 version*

This workflow uses the capability of mutect2 to call several samples at the same time and improve recall and accuracy through a joint model.
        Most of these tools are still in a beta state and not intended for main production (as of 4.1.4.0)
        There are also som major tweaks we have to do for runtime, as the amount of data might overwhelm the tools otherwise.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.dawson.workflows.mutectjointsomaticworkflow import Mutect2JointSomaticWorkflow

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "mutect2jointsomaticworkflow_step",
           Mutect2JointSomaticWorkflow(
               normalBams=None,
               tumorBams=None,
               normalName=None,
               biallelicSites=None,
               reference=None,
               panelOfNormals=None,
               germlineResource=None,
           )
       )
       wf.output("out", source=mutect2jointsomaticworkflow_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Mutect2JointSomaticWorkflow:

.. code-block:: bash

   # user inputs
   janis inputs Mutect2JointSomaticWorkflow > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       biallelicSites: biallelicSites.vcf.gz
       germlineResource: germlineResource.vcf.gz
       normalBams:
       - normalBams_0.cram
       - normalBams_1.cram
       normalName: <value>
       panelOfNormals: panelOfNormals.vcf.gz
       reference: reference.fasta
       tumorBams:
       - tumorBams_0.cram
       - tumorBams_1.cram




5. Run Mutect2JointSomaticWorkflow with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Mutect2JointSomaticWorkflow





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``Mutect2JointSomaticWorkflow``
:URL: *No URL to the documentation was provided*
:Versions: 0.1
:Authors: 
:Citations: 
:Created: 2019-10-30
:Updated: 2019-10-30



Outputs
-----------

======  ====================  ===============
name    type                  documentation
======  ====================  ===============
out     CompressedIndexedVCF
======  ====================  ===============


Embedded Tools
***************

================================  ==========================================
Create genomic call regions       ``CreateCallRegions/v0.1.0``
GatkMutect2                       ``Gatk4Mutect2_cram/4.1.4.0``
BCFTools: Concat                  ``bcftoolsConcat/v1.9``
BCFTools: Index                   ``bcftoolsIndex/v1.9``
GATK4: LearnReadOrientationModel  ``Gatk4LearnReadOrientationModel/4.1.4.0``
GATK4: MergeMutectStats           ``Gatk4MergeMutectStats/4.1.2.0``
GATK4: GetPileupSummaries         ``Gatk4GetPileupSummaries_cram/4.1.4.0``
GATK4: CalculateContamination     ``Gatk4CalculateContamination/4.1.4.0``
GATK4: GetFilterMutectCalls       ``Gatk4FilterMutectCalls/4.1.4.0``
BCFTools: Normalize               ``bcftoolsNorm/v1.9``
================================  ==========================================



Additional configuration (inputs)
---------------------------------

==========================  ====================  ===============
name                        type                  documentation
==========================  ====================  ===============
normalBams                  Array<CramPair>
tumorBams                   Array<CramPair>
normalName                  String
biallelicSites              CompressedIndexedVCF
reference                   FastaWithIndexes
panelOfNormals              CompressedIndexedVCF
germlineResource            CompressedIndexedVCF
regionSize                  Optional<Integer>
createCallRegions_equalize  Optional<Boolean>
==========================  ====================  ===============


