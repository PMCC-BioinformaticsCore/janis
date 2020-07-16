:orphan:

BWA-Index
====================

*1 contributor · 1 version*

bwa - Burrows-Wheeler Alignment Tool
Index database sequences in the FASTA format.

Warning: `-a bwtsw' does not work for short genomes, while `-a is' and
         `-a div' do not work not for long genomes.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.bwa.index.versions import BwaIndex_0_7_15

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bwaindex_step",
           BwaIndex_0_7_15(
               reference=None,
           )
       )
       wf.output("out", source=bwaindex_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bwaIndex:

.. code-block:: bash

   # user inputs
   janis inputs bwaIndex > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run bwaIndex with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bwaIndex





Information
------------


:ID: ``bwaIndex``
:URL: `http://bio-bwa.sourceforge.net/bwa.shtml#3 <http://bio-bwa.sourceforge.net/bwa.shtml#3>`_
:Versions: v0.7.15
:Container: biocontainers/bwa:v0.7.15_cv3
:Authors: Michael Franklin
:Citations: The BWA-MEM algorithm has not been published yet.
:Created: 2020-02-14
:Updated: 2020-02-14



Outputs
-----------

======  ========  ===============
name    type      documentation
======  ========  ===============
out     FastaBwa
======  ========  ===============



Additional configuration (inputs)
---------------------------------

=========  =================  ========  ==========  =======================================================================
name       type               prefix      position  documentation
=========  =================  ========  ==========  =======================================================================
reference  Fasta                                 1
blockSize  Optional<Integer>  -b                    block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]
algorithm  Optional<String>   -a                    BWT construction algorithm: bwtsw, is or rb2 [auto]
=========  =================  ========  ==========  =======================================================================
