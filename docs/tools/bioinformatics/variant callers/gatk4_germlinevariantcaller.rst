:orphan:

GATK4 Germline Variant Caller
===========================================================

*1 contributor · 2 versions*

This is a VariantCaller based on the GATK Best Practice pipelines. It uses the GATK4 toolkit, specifically 4.0.12.0.

        It has the following steps:

        1. BaseRecalibrator
        2. ApplyBQSR
        3. HaplotypeCaller
        4. SplitMultiAllele


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.variantcallers.gatk.gatkgermline_variants_4_1_3 import GatkGermlineVariantCaller_4_1_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4_germlinevariantcaller_step",
           GatkGermlineVariantCaller_4_1_3(
               bam=None,
               reference=None,
               snps_dbsnp=None,
               snps_1000gp=None,
               known_indels=None,
               mills_indels=None,
           )
       )
       wf.output("out", source=gatk4_germlinevariantcaller_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for GATK4_GermlineVariantCaller:

.. code-block:: bash

   # user inputs
   janis inputs GATK4_GermlineVariantCaller > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       known_indels: known_indels.vcf.gz
       mills_indels: mills_indels.vcf.gz
       reference: reference.fasta
       snps_1000gp: snps_1000gp.vcf.gz
       snps_dbsnp: snps_dbsnp.vcf.gz




5. Run GATK4_GermlineVariantCaller with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       GATK4_GermlineVariantCaller





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``GATK4_GermlineVariantCaller``
:URL: *No URL to the documentation was provided*
:Versions: 4.0.12.0, 4.1.3.0
:Authors: Michael Franklin
:Citations: 
:Created: 2019-09-01
:Updated: 2019-09-13



Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============


Embedded Tools
***************

=============================================  =================================
GATK4: SplitReads                              ``Gatk4SplitReads/4.1.3.0``
GATK4: Base Recalibrator                       ``Gatk4BaseRecalibrator/4.1.3.0``
GATK4: Apply base quality score recalibration  ``Gatk4ApplyBQSR/4.1.3.0``
GATK4: Haplotype Caller                        ``Gatk4HaplotypeCaller/4.1.3.0``
Split Multiple Alleles                         ``SplitMultiAllele/v0.5772``
=============================================  =================================



Additional configuration (inputs)
---------------------------------

============  ====================  ===================================================================================================================================================
name          type                  documentation
============  ====================  ===================================================================================================================================================
bam           IndexedBam
reference     FastaWithIndexes
snps_dbsnp    CompressedIndexedVCF
snps_1000gp   CompressedIndexedVCF
known_indels  CompressedIndexedVCF
mills_indels  CompressedIndexedVCF
intervals     Optional<bed>         This optional interval supports processing by regions. If this input resolves to null, then GATK will process the whole genome per each tool's spec
============  ====================  ===================================================================================================================================================


