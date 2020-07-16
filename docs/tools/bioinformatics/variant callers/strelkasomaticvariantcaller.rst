:orphan:

Strelka Somatic Variant Caller
============================================================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.variantcallers.illuminasomatic_strelka import IlluminaSomaticVariantCaller

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "strelkasomaticvariantcaller_step",
           IlluminaSomaticVariantCaller(
               normal_bam=None,
               tumor_bam=None,
               reference=None,
           )
       )
       wf.output("sv", source=strelkasomaticvariantcaller_step.sv)
       wf.output("variants", source=strelkasomaticvariantcaller_step.variants)
       wf.output("out", source=strelkasomaticvariantcaller_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for strelkaSomaticVariantCaller:

.. code-block:: bash

   # user inputs
   janis inputs strelkaSomaticVariantCaller > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       normal_bam: normal_bam.bam
       reference: reference.fasta
       tumor_bam: tumor_bam.bam




5. Run strelkaSomaticVariantCaller with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       strelkaSomaticVariantCaller





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``strelkaSomaticVariantCaller``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.1
:Authors: 
:Citations: 
:Created: None
:Updated: None



Outputs
-----------

========  ====================  ===============
name      type                  documentation
========  ====================  ===============
sv        CompressedIndexedVCF
variants  CompressedIndexedVCF
out       VCF
========  ====================  ===============


Embedded Tools
***************

==========================  ==================================
Manta                       ``manta/1.5.0``
Strelka (Somatic)           ``strelka_somatic/2.9.10``
Concat Strelka Somatic Vcf  ``ConcatStrelkaSomaticVcf/0.1.16``
Tabix                       ``tabix/1.2.1``
UncompressArchive           ``UncompressArchive/v1.0.0``
Split Multiple Alleles      ``SplitMultiAllele/v0.5772``
VcfTools                    ``VcfTools/0.1.16``
==========================  ==================================



Additional configuration (inputs)
---------------------------------

=============================  ==================  =================================================================================================================================================================================================================================================================
name                           type                documentation
=============================  ==================  =================================================================================================================================================================================================================================================================
normal_bam                     IndexedBam
tumor_bam                      IndexedBam
reference                      FastaWithIndexes
intervals                      Optional<BedTABIX>
is_exome                       Optional<Boolean>
filterpass_removeFileteredAll  Optional<Boolean>   Removes all sites with a FILTER flag other than PASS.
filterpass_recode              Optional<Boolean>
filterpass_recodeINFOAll       Optional<Boolean>   These options can be used with the above recode options to define an INFO key name to keep in the output  file.  This  option can be used multiple times to keep more of the INFO fields. The second option is used to keep all INFO values in the original file.
=============================  ==================  =================================================================================================================================================================================================================================================================


