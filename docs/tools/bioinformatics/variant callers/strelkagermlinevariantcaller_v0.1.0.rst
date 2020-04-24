:orphan:

Strelka Germline Variant Caller
==============================================================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.variantcallers.illuminagermline_strelka import IlluminaGermlineVariantCaller

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "strelkagermlinevariantcaller_step",
           IlluminaGermlineVariantCaller(
               bam=None,
               reference=None,
           )
       )
       wf.output("diploid", source=strelkagermlinevariantcaller_step.diploid)
   wf.output("variants", source=strelkagermlinevariantcaller_step.variants)
   wf.output("out", source=strelkagermlinevariantcaller_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for strelkaGermlineVariantCaller:

.. code-block:: bash

   # user inputs
   janis inputs strelkaGermlineVariantCaller > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       reference: reference.fasta




5. Run strelkaGermlineVariantCaller with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       strelkaGermlineVariantCaller





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``strelkaGermlineVariantCaller``
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
Strelka (Germline)      ``strelka_germline/2.9.10``
BCFTools: View          ``bcftoolsview/v1.5``
Split Multiple Alleles  ``SplitMultiAllele/v0.5772``
======================  ============================



Additional configuration (inputs)
---------------------------------

====================  =======================  =======================================================================
name                  type                     documentation
====================  =======================  =======================================================================
bam                   IndexedBam
reference             FastaWithIndexes
intervals             Optional<BedTABIX>
is_exome              Optional<Boolean>
bcfview_applyFilters  Optional<Array<String>>  (-f) require at least one of the listed FILTER strings (e.g. 'PASS,.'')
====================  =======================  =======================================================================


