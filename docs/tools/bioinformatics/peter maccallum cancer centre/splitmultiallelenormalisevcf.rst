:orphan:

Split Multiple Alleles and Normalise Vcf
=======================================================================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.common.splitmultiallele_normalistvcf import SplitMultiAlleleNormaliseVcf

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "splitmultiallelenormalisevcf_step",
           SplitMultiAlleleNormaliseVcf(
               reference=None,
           )
       )
       wf.output("out", source=splitmultiallelenormalisevcf_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SplitMultiAlleleNormaliseVcf:

.. code-block:: bash

   # user inputs
   janis inputs SplitMultiAlleleNormaliseVcf > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run SplitMultiAlleleNormaliseVcf with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SplitMultiAlleleNormaliseVcf





Information
------------


:ID: ``SplitMultiAlleleNormaliseVcf``
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

==================  ==============================  ========  ==========  ===============
name                type                            prefix      position  documentation
==================  ==============================  ========  ==========  ===============
reference           FastaWithIndexes                -r                 4
vcf                 Optional<VCF>                                      1
compressedTabixVcf  Optional<CompressedIndexedVCF>                     1
compressedVcf       Optional<CompressedVCF>                            1
outputFilename      Optional<Filename>              -o                 6
==================  ==============================  ========  ==========  ===============
