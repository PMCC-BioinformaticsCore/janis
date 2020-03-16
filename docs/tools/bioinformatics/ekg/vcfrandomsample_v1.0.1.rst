:orphan:

VcfLib: Vcf Random Sampling
=============================================

*1 contributor · 1 version*

usage: vcfrandomsample [options] [<vcf file>]

options:
	-r, --rate RATE 	base sampling probability per locus
	-s, --scale-by KEY\scale sampling likelihood by this Float info field
	-p, --random-seed N	use this random seed

Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcfrandomsample.versions import VcfRandomSample_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcfrandomsample_step",
           vcfrandomsample(
               vcf=None,
               rate=None,
               seed=None,
           )
       )
       wf.output("out", source=vcfrandomsample_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcfrandomsample:

.. code-block:: bash

   # user inputs
   janis inputs vcfrandomsample > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       rate: 0.0
       seed: 0
       vcf: vcf.vcf.gz




5. Run vcfrandomsample with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcfrandomsample





Information
------------


:ID: ``vcfrandomsample``
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

=======  ================  ========  ==========  ==================================================
name     type              prefix      position  documentation
=======  ================  ========  ==========  ==================================================
vcf      CompressedVCF                        3
rate     Float             -t                    base sampling probability per locus
seed     Integer           -p                    use this random seed
scaleBy  Optional<String>  -s                    scale sampling likelihood by this Float info field
=======  ================  ========  ==========  ==================================================
