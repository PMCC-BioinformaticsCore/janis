:orphan:

VcfLib: VcfUniqAlleles
=======================================

*1 contributor · 1 version*

usage: vcffuniq [file]
For each record, remove any duplicate alternate alleles that may have resulted from merging separate VCF files.

Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcfuniqalleles.versions import VcfUniqAlleles_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcfuniqalleles_step",
           vcfuniqalleles(
               vcf=None,
           )
       )
       wf.output("out", source=vcfuniqalleles_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcfuniqalleles:

.. code-block:: bash

   # user inputs
   janis inputs vcfuniqalleles > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf.gz




5. Run vcfuniqalleles with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcfuniqalleles





Information
------------


:ID: ``vcfuniqalleles``
:URL: `https://github.com/vcflib/vcflib <https://github.com/vcflib/vcflib>`_
:Versions: v1.0.1
:Container: shollizeck/vcflib:1.0.1
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2019-10-18
:Updated: 2019-10-18



Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     stdout<VCF>  VCF output
======  ===========  ===============



Additional configuration (inputs)
---------------------------------

======  =============  ========  ==========  ===============
name    type           prefix      position  documentation
======  =============  ========  ==========  ===============
vcf     CompressedVCF                     3
======  =============  ========  ==========  ===============
