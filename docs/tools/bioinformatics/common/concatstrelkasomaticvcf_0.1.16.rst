:orphan:

Concat Strelka Somatic Vcf
====================================================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.common.concat_strelkasomaticvcf import ConcatStrelkaSomaticVcf

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "concatstrelkasomaticvcf_step",
           ConcatStrelkaSomaticVcf(
               headerVcfs=None,
               contentVcfs=None,
           )
       )
       wf.output("out", source=concatstrelkasomaticvcf_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for ConcatStrelkaSomaticVcf:

.. code-block:: bash

   # user inputs
   janis inputs ConcatStrelkaSomaticVcf > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       contentVcfs:
       - contentVcfs_0.vcf.gz
       - contentVcfs_1.vcf.gz
       headerVcfs:
       - headerVcfs_0.vcf.gz
       - headerVcfs_1.vcf.gz




5. Run ConcatStrelkaSomaticVcf with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       ConcatStrelkaSomaticVcf





Information
------------


:ID: ``ConcatStrelkaSomaticVcf``
:URL: *No URL to the documentation was provided*
:Versions: 0.1.16
:Container: biocontainers/vcftools:v0.1.16-1-deb_cv1
:Authors: 
:Citations: None
:Created: None
:Updated: None



Outputs
-----------

======  =============  ===============
name    type           documentation
======  =============  ===============
out     CompressedVCF
======  =============  ===============



Additional configuration (inputs)
---------------------------------

==============  ===========================  ========  ==========  ===============
name            type                         prefix      position  documentation
==============  ===========================  ========  ==========  ===============
headerVcfs      Array<CompressedIndexedVCF>                     1
contentVcfs     Array<CompressedIndexedVCF>                     4
outputFilename  Optional<Filename>           >                  6
==============  ===========================  ========  ==========  ===============
