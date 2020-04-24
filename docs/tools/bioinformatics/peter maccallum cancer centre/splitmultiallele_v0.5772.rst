:orphan:

Split Multiple Alleles
=========================================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.common.splitmultiallele import SplitMultiAllele

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "splitmultiallele_step",
           SplitMultiAllele(
               vcf=None,
               reference=None,
           )
       )
       wf.output("out", source=splitmultiallele_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SplitMultiAllele:

.. code-block:: bash

   # user inputs
   janis inputs SplitMultiAllele > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta
       vcf: vcf.vcf.gz




5. Run SplitMultiAllele with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SplitMultiAllele





Information
------------


:ID: ``SplitMultiAllele``
:URL: *No URL to the documentation was provided*
:Versions: v0.5772
:Container: heuermh/vt
:Authors: 
:Citations: None
:Created: None
:Updated: None



Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============



Additional configuration (inputs)
---------------------------------

==============  ==================  ========  ==========  ===============
name            type                prefix      position  documentation
==============  ==================  ========  ==========  ===============
vcf             CompressedVCF                          3
reference       FastaWithIndexes    -r                 8
outputFilename  Optional<Filename>  >                 10
==============  ==================  ========  ==========  ===============
