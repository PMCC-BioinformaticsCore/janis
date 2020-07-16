:orphan:

VcfLib: VcfStreamSort
=====================================

*1 contributor · 1 version*

usage: vcfallelicprimitives [options] [file]

options:
	-m, --use-mnps	Retain MNPs as separate events (default: false)
	-t, --tag-parsed FLAG	Tag records which are split apart of a complex allele with this flag


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcfstreamsort.versions import VcfStreamSort_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcfstreamsort_step",
           VcfStreamSort_1_0_1(
               vcf=None,
           )
       )
       wf.output("out", source=vcfstreamsort_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcfstreamsort:

.. code-block:: bash

   # user inputs
   janis inputs vcfstreamsort > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf




5. Run vcfstreamsort with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcfstreamsort





Information
------------


:ID: ``vcfstreamsort``
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

============  =================  ========  ==========  =======================================
name          type               prefix      position  documentation
============  =================  ========  ==========  =======================================
vcf           VCF                                   3
inMemoryFlag  Optional<Boolean>  -a                    load all sites and then sort in memory
windowSize    Optional<Integer>  -w                    number of sites to sort (default 10000)
============  =================  ========  ==========  =======================================
