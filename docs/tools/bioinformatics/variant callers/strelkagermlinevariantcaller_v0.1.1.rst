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
       wf.output("sv", source=strelkagermlinevariantcaller_step.sv)
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

======================  ============================
Manta                   ``manta/1.5.0``
Strelka (Germline)      ``strelka_germline/2.9.10``
UncompressArchive       ``UncompressArchive/v1.0.0``
Split Multiple Alleles  ``SplitMultiAllele/v0.5772``
VcfTools                ``VcfTools/0.1.16``
======================  ============================



Additional configuration (inputs)
---------------------------------

=============================  ==================  =================================================================================================================================================================================================================================================================
name                           type                documentation
=============================  ==================  =================================================================================================================================================================================================================================================================
bam                            IndexedBam
reference                      FastaWithIndexes
intervals                      Optional<BedTABIX>
is_exome                       Optional<Boolean>
filterpass_removeFileteredAll  Optional<Boolean>   Removes all sites with a FILTER flag other than PASS.
filterpass_recode              Optional<Boolean>
filterpass_recodeINFOAll       Optional<Boolean>   These options can be used with the above recode options to define an INFO key name to keep in the output  file.  This  option can be used multiple times to keep more of the INFO fields. The second option is used to keep all INFO values in the original file.
=============================  ==================  =================================================================================================================================================================================================================================================================


