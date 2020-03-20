:orphan:

VcfLib: Vcf ROC generator
==================================

*1 contributor · 1 version*

usage: vcfroc [options] [<vcf file>]

options:
	-t, --truth-vcf FILE	use this VCF as ground truth for ROC generation
	-w, --window-size N       compare records up to this many bp away (default 30)
	-r, --reference FILE	FASTA reference file

Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcfroc.versions import VcfRoc_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcfroc_step",
           VcfRoc_1_0_1(
               vcf=None,
               truth=None,
               reference=None,
           )
       )
       wf.output("out", source=vcfroc_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcfroc:

.. code-block:: bash

   # user inputs
   janis inputs vcfroc > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta
       truth: truth.vcf.gz
       vcf: vcf.vcf.gz




5. Run vcfroc with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcfroc





Information
------------


:ID: ``vcfroc``
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

==========  =================  ========  ==========  ====================================================
name        type               prefix      position  documentation
==========  =================  ========  ==========  ====================================================
vcf         CompressedVCF                         3
truth       CompressedVCF      -t                    use this VCF as ground truth for ROC generation
reference   FastaWithIndexes   -r                    FASTA reference file
windowSize  Optional<Integer>  -w                    compare records up to this many bp away (default 30)
==========  =================  ========  ==========  ====================================================
