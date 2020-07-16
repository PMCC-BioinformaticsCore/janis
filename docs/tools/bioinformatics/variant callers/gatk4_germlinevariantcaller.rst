:orphan:

GATK4 Germline Variant Caller
===========================================================

*2 contributors · 2 versions*

This is a VariantCaller based on the GATK Best Practice pipelines. It uses the GATK4 toolkit, specifically 4.1.3.

        It has the following steps:

        1. Split Bam based on intervals (bed)
        2. HaplotypeCaller
        3. SplitMultiAllele


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
           )
       )
       wf.output("variants", source=gatk4_germlinevariantcaller_step.variants)
       wf.output("out_bam", source=gatk4_germlinevariantcaller_step.out_bam)
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
       reference: reference.fasta
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
:Authors: Michael Franklin, Jiaan
:Citations: 
:Created: 2019-09-01
:Updated: 2019-09-13



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

=======================  ================================
GATK4: SplitReads        ``Gatk4SplitReads/4.1.3.0``
GATK4: Haplotype Caller  ``Gatk4HaplotypeCaller/4.1.3.0``
UncompressArchive        ``UncompressArchive/v1.0.0``
Split Multiple Alleles   ``SplitMultiAllele/v0.5772``
=======================  ================================



Additional configuration (inputs)
---------------------------------

======================================  ====================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================
name                                    type                  documentation
======================================  ====================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================
bam                                     IndexedBam
reference                               FastaWithIndexes
snps_dbsnp                              CompressedIndexedVCF
intervals                               Optional<bed>         This optional interval supports processing by regions. If this input resolves to null, then GATK will process the whole genome per each tool's spec
haplotype_caller_pairHmmImplementation  Optional<String>      The PairHMM implementation to use for genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime. The --pair-hmm-implementation argument is an enumerated type (Implementation), which can have one of the following values: EXACT;ORIGINAL;LOGLESS_CACHING;AVX_LOGLESS_CACHING;AVX_LOGLESS_CACHING_OMP;EXPERIMENTAL_FPGA_LOGLESS_CACHING;FASTEST_AVAILABLE. Implementation:  FASTEST_AVAILABLE
======================================  ====================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================


