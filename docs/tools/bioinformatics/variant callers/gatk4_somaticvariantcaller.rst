:orphan:

GATK4 Somatic Variant Caller
=========================================================

*2 contributors · 2 versions*

This is a VariantCaller based on the GATK Best Practice pipelines. It uses the GATK4 toolkit, specifically 4.0.12.0. Takes GATK Base Recalibrated Bam as input

        It has the following steps:

        1. Mutect2
        2. LearnOrientationModel
        3. GetPileUpSummaries
        4. CalculateContamination
        5. FilterMutectCalls
        6. Split and normliase vcf
        7. Filter PASS variants


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.variantcallers.gatk.gatksomatic_variants_4_1_3 import GatkSomaticVariantCaller_4_1_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4_somaticvariantcaller_step",
           GatkSomaticVariantCaller_4_1_3(
               normal_bam=None,
               tumor_bam=None,
               reference=None,
               gnomad=None,
           )
       )
       wf.output("variants", source=gatk4_somaticvariantcaller_step.variants)
       wf.output("out_bam", source=gatk4_somaticvariantcaller_step.out_bam)
       wf.output("out", source=gatk4_somaticvariantcaller_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for GATK4_SomaticVariantCaller:

.. code-block:: bash

   # user inputs
   janis inputs GATK4_SomaticVariantCaller > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       gnomad: gnomad.vcf.gz
       normal_bam: normal_bam.bam
       reference: reference.fasta
       tumor_bam: tumor_bam.bam




5. Run GATK4_SomaticVariantCaller with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       GATK4_SomaticVariantCaller





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``GATK4_SomaticVariantCaller``
:URL: *No URL to the documentation was provided*
:Versions: 4.0.12.0, 4.1.3.0
:Authors: Michael Franklin, Jiaan Yu
:Citations: 
:Created: 2019-02-01
:Updated: 2020-06-15



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

================================  ==========================================
                                  ``split_bam_subpipeline/None``
GatkMutect2                       ``Gatk4Mutect2/4.1.3.0``
GATK4: LearnReadOrientationModel  ``Gatk4LearnReadOrientationModel/4.1.4.0``
GATK4: GetPileupSummaries         ``Gatk4GetPileupSummaries/4.1.4.0``
GATK4: CalculateContamination     ``Gatk4CalculateContamination/4.1.4.0``
GATK4: GetFilterMutectCalls       ``Gatk4FilterMutectCalls/4.1.3.0``
UncompressArchive                 ``UncompressArchive/v1.0.0``
Split Multiple Alleles            ``SplitMultiAllele/v0.5772``
VcfTools                          ``VcfTools/0.1.16``
================================  ==========================================



Additional configuration (inputs)
---------------------------------

=============================  ==============================  =================================================================================================================================================================================================================================================================
name                           type                            documentation
=============================  ==============================  =================================================================================================================================================================================================================================================================
normal_bam                     IndexedBam
tumor_bam                      IndexedBam
reference                      FastaWithIndexes
gnomad                         CompressedIndexedVCF
normal_name                    Optional<String>
intervals                      Optional<bed>                   This optional intervals file supports processing by regions. If this file resolves to null, then GATK will process the whole genome per each tool's spec
panel_of_normals               Optional<CompressedIndexedVCF>
filterpass_removeFileteredAll  Optional<Boolean>               Removes all sites with a FILTER flag other than PASS.
filterpass_recode              Optional<Boolean>
filterpass_recodeINFOAll       Optional<Boolean>               These options can be used with the above recode options to define an INFO key name to keep in the output  file.  This  option can be used multiple times to keep more of the INFO fields. The second option is used to keep all INFO values in the original file.
=============================  ==============================  =================================================================================================================================================================================================================================================================


