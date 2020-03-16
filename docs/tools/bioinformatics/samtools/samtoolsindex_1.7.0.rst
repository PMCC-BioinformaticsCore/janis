:orphan:

SamTools: Index
===============================

*1 contributor · 2 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.samtools.index.versions import SamToolsIndex_1_7

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "samtoolsindex_step",
           SamToolsIndex(
               bam=None,
           )
       )
       wf.output("out", source=samtoolsindex_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SamToolsIndex:

.. code-block:: bash

   # user inputs
   janis inputs SamToolsIndex > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam




5. Run SamToolsIndex with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SamToolsIndex





Information
------------


:ID: ``SamToolsIndex``
:URL: `http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS <http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS>`_
:Versions: 1.9.0, 1.7.0
:Container: biocontainers/samtools:v1.7.0_cv3
:Authors: Michael Franklin
:Citations: None
:Created: 2019-12-17
:Updated: 2019-12-17



Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedBam
======  ==========  ===============



Additional configuration (inputs)
---------------------------------

=======  =================  ========  ==========  ===============
name     type               prefix      position  documentation
=======  =================  ========  ==========  ===============
bam      BAM                                  10
threads  Optional<Integer>  -@                10
=======  =================  ========  ==========  ===============
