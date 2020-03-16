:orphan:

GATK4: CreateSequenceDictionary
===============================================================

*1 contributor · 3 versions*

Creates a sequence dictionary for a reference sequence.  This tool creates a sequence dictionary file (with ".dict"
extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools.
The output file contains a header but no SAMRecords, and the header contains only sequence records.

The reference sequence can be gzipped (both .fasta and .fasta.gz are supported).

Usage example:

    java -jar picard.jar CreateSequenceDictionary \
        R=reference.fasta \
        O=reference.dict

Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.createsequencedictionary.versions import Gatk4CreateSequenceDictionary_4_1_4

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4createsequencedictionary_step",
           Gatk4CreateSequenceDictionary(
               reference=None,
           )
       )
       wf.output("out", source=gatk4createsequencedictionary_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4CreateSequenceDictionary:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4CreateSequenceDictionary > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run Gatk4CreateSequenceDictionary with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4CreateSequenceDictionary





Information
------------


:ID: ``Gatk4CreateSequenceDictionary``
:URL: `https://gatk.broadinstitute.org/hc/en-us/articles/360036509572-CreateSequenceDictionary-Picard- <https://gatk.broadinstitute.org/hc/en-us/articles/360036509572-CreateSequenceDictionary-Picard->`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.4.0
:Authors: Michael Franklin
:Citations: TBD
:Created: 2020-02-14
:Updated: 2020-02-14



Outputs
-----------

======  ========  ======================================
name    type      documentation
======  ========  ======================================
out     FastDict  Output reference with ^.dict reference
======  ========  ======================================



Additional configuration (inputs)
---------------------------------

=========  ======  ===========  ==========  =================================================
name       type    prefix       position    documentation
=========  ======  ===========  ==========  =================================================
reference  Fasta   --REFERENCE              (-R) Input reference fasta or fasta.gz  Required.
=========  ======  ===========  ==========  =================================================
