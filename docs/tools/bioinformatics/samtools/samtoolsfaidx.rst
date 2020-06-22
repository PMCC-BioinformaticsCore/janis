:orphan:

SamTools: faidx
===============================

*1 contributor · 2 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.samtools.faidx.versions import SamToolsFaidx_1_9

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "samtoolsfaidx_step",
           SamToolsFaidx_1_9(
               reference=None,
           )
       )
       wf.output("out", source=samtoolsfaidx_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SamToolsFaidx:

.. code-block:: bash

   # user inputs
   janis inputs SamToolsFaidx > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run SamToolsFaidx with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SamToolsFaidx





Information
------------


:ID: ``SamToolsFaidx``
:URL: `http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS <http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS>`_
:Versions: 1.9.0, 1.7.0
:Container: quay.io/biocontainers/samtools:1.9--h8571acd_11
:Authors: Michael Franklin
:Citations: None
:Created: 2020-02-14
:Updated: 2020-02-14



Outputs
-----------

======  ========  ===============
name    type      documentation
======  ========  ===============
out     FastaFai
======  ========  ===============



Additional configuration (inputs)
---------------------------------

=========  ======  ========  ==========  ===============
name       type    prefix      position  documentation
=========  ======  ========  ==========  ===============
reference  Fasta                      1
=========  ======  ========  ==========  ===============
