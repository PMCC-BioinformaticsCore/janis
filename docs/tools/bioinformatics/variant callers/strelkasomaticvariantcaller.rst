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
       wf.output("diploid", source=strelkasomaticvariantcaller_step.diploid)
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
:Versions: v0.1.0
:Authors: 
:Citations: 
:Created: None
:Updated: None



Outputs
-----------

========  ====================  ===============
name      type                  documentation
========  ====================  ===============
diploid   CompressedIndexedVCF
variants  CompressedIndexedVCF
out       VCF
========  ====================  ===============


Embedded Tools
***************

======================  ============================
Manta                   ``manta/1.5.0``
Strelka (Somatic)       ``strelka_somatic/2.9.10``
BCFTools: View          ``bcftoolsview/v1.5``
Split Multiple Alleles  ``SplitMultiAllele/v0.5772``
======================  ============================



Additional configuration (inputs)
---------------------------------

=====================  =======================  =======================================================================
name                   type                     documentation
=====================  =======================  =======================================================================
normal_bam             IndexedBam
tumor_bam              IndexedBam
reference              FastaWithIndexes
intervals              Optional<BedTABIX>
is_exome               Optional<Boolean>
bcf_view_applyFilters  Optional<Array<String>>  (-f) require at least one of the listed FILTER strings (e.g. 'PASS,.'')
=====================  =======================  =======================================================================


