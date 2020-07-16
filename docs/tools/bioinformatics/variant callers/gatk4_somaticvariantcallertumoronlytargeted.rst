:orphan:

GATK4 Somatic Variant Caller for Tumour Only Samples with Targeted BED
====================================================================================================================

*2 contributors · 1 version*

This is a VariantCaller based on the GATK Best Practice pipelines. It uses the GATK4 toolkit, specifically 4.1.2.

        It has the following steps:

        1. Mutect2 (output: vcf, bam, f1r2.tar.gz)
        2. LearnOrientationModel
        3. GetPileupSummaries
        4. CalculateContamination
        5. FilterMutect2Calls
        6. SplitNormaliseVcf


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.variantcallers.gatk.gatksomatic_variants_single import GatkSomaticVariantCallerTumorOnlyTargeted

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4_somaticvariantcallertumoronlytargeted_step",
           GatkSomaticVariantCallerTumorOnlyTargeted(
               bam=None,
               reference=None,
               gnomad=None,
           )
       )
       wf.output("variants", source=gatk4_somaticvariantcallertumoronlytargeted_step.variants)
       wf.output("out_bam", source=gatk4_somaticvariantcallertumoronlytargeted_step.out_bam)
       wf.output("out", source=gatk4_somaticvariantcallertumoronlytargeted_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for GATK4_SomaticVariantCallerTumorOnlyTargeted:

.. code-block:: bash

   # user inputs
   janis inputs GATK4_SomaticVariantCallerTumorOnlyTargeted > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       gnomad: gnomad.vcf.gz
       reference: reference.fasta




5. Run GATK4_SomaticVariantCallerTumorOnlyTargeted with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       GATK4_SomaticVariantCallerTumorOnlyTargeted





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``GATK4_SomaticVariantCallerTumorOnlyTargeted``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.1
:Authors: Michael Franklin, Jiaan Yu
:Citations: 
:Created: 2020-06-04
:Updated: 2020-06-04



Outputs
-----------

========  ====================  ===============
name      type                  documentation
========  ====================  ===============
variants  CompressedIndexedVCF
out_bam   IndexedBam
out       VCF
========  ====================  ===============


Embedded Tools
***************

========================================  ==========================================
GatkMutect2                               ``Gatk4Mutect2/4.1.2.0``
GATK4: LearnReadOrientationModel          ``Gatk4LearnReadOrientationModel/4.1.2.0``
GATK4: GetPileupSummaries                 ``Gatk4GetPileupSummaries/4.1.2.0``
GATK4: CalculateContamination             ``Gatk4CalculateContamination/4.1.2.0``
GATK4: GetFilterMutectCalls               ``Gatk4FilterMutectCalls/4.1.2.0``
Split Multiple Alleles and Normalise Vcf  ``SplitMultiAlleleNormaliseVcf/v0.5772``
========================================  ==========================================



Additional configuration (inputs)
---------------------------------

================  ==============================  ===================================================================================================================================================
name              type                            documentation
================  ==============================  ===================================================================================================================================================
bam               IndexedBam
reference         FastaWithIndexes
gnomad            CompressedIndexedVCF
intervals         Optional<bed>                   This optional interval supports processing by regions. If this input resolves to null, then GATK will process the whole genome per each tool's spec
panel_of_normals  Optional<CompressedIndexedVCF>
================  ==============================  ===================================================================================================================================================


