:orphan:

VcfLib: VcfCombine
===============================

*1 contributor · 1 version*

usage: vcfcombine [vcf file] [vcf file] ...

 options:
-h --help	This text.
-r --region REGION	A region specifier of the form chrN:x-y to bound the merge

Combines VCF files positionally, combining samples when sites and alleles are identical. Any number of VCF files may be combined. The INFO field and other columns are taken from one of the files which are combined when records in multiple files match. Alleles must have identical ordering to be combined into one record. If they do not, multiple records will be emitted


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcfcombine.versions import VcfCombine_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcfcombine_step",
           VcfCombine_1_0_1(
               vcf=None,
           )
       )
       wf.output("out", source=vcfcombine_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcfcombine:

.. code-block:: bash

   # user inputs
   janis inputs vcfcombine > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf:
       - vcf_0.vcf
       - vcf_1.vcf




5. Run vcfcombine with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcfcombine





Information
------------


:ID: ``vcfcombine``
:URL: `https://github.com/vcflib/vcflib <https://github.com/vcflib/vcflib>`_
:Versions: v1.0.1
:Container: shollizeck/vcflib:1.0.1
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2019-10-18
:Updated: 2019-10-18



Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     stdout<File>  VCF output
======  ============  ===============



Additional configuration (inputs)
---------------------------------

======  ================  ========  ==========  ==========================================================
name    type              prefix      position  documentation
======  ================  ========  ==========  ==========================================================
vcf     Array<VCF>                           2
region  Optional<String>  -r                 1  A region specifier of the form chrN:x-y to bound the merge
======  ================  ========  ==========  ==========================================================
