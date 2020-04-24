:orphan:

VcfLib: VcfUniq
=========================

*1 contributor · 1 version*

usage: vcffuniq [file]
Like GNU uniq, but for VCF records. Remove records which have the same positon, ref, and alt as the previous record.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcfuniq.versions import VcfUniq_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcfuniq_step",
           VcfUniq_1_0_1(
               vcf=None,
           )
       )
       wf.output("out", source=vcfuniq_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcfuniq:

.. code-block:: bash

   # user inputs
   janis inputs vcfuniq > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf.gz




5. Run vcfuniq with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcfuniq





Information
------------


:ID: ``vcfuniq``
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
