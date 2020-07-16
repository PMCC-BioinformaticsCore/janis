:orphan:

VcfLib: Vcf Length
==============================

*1 contributor · 1 version*

Adds the length of the variant record (in [-/+]) relative to the reference allele to each VCF record.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcflength.versions import VcfLength_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcflength_step",
           VcfLength_1_0_1(
               vcf=None,
           )
       )
       wf.output("out", source=vcflength_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcflength:

.. code-block:: bash

   # user inputs
   janis inputs vcflength > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf




5. Run vcflength with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcflength





Information
------------


:ID: ``vcflength``
:URL: `https://github.com/vcflib/vcflib <https://github.com/vcflib/vcflib>`_
:Versions: v1.0.1
:Container: shollizeck/vcflib:1.0.1
:Authors: Michael Franklin
:Citations: None
:Created: 2020-06-04
:Updated: 2020-06-04



Outputs
-----------

======  ===========  ==========================================================
name    type         documentation
======  ===========  ==========================================================
out     stdout<VCF>  VCF with length of the variant record added to each record
======  ===========  ==========================================================



Additional configuration (inputs)
---------------------------------

======  ======  ========  ==========  ========================================================================
name    type    prefix      position  documentation
======  ======  ========  ==========  ========================================================================
vcf     VCF                        1  VCF to add length of variant record relative to the reference allele to.
======  ======  ========  ==========  ========================================================================
