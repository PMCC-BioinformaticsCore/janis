:orphan:

Index Fasta reference
==================================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.common.indexfasta import IndexFasta

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "indexfasta_step",
           IndexFasta(
               reference=None,
           )
       )
       wf.output("bwa", source=indexfasta_step.bwa)
       wf.output("samtools", source=indexfasta_step.samtools)
       wf.output("dict", source=indexfasta_step.dict)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for IndexFasta:

.. code-block:: bash

   # user inputs
   janis inputs IndexFasta > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run IndexFasta with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       IndexFasta





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``IndexFasta``
:URL: *No URL to the documentation was provided*
:Versions: 1.0.0
:Authors: 
:Citations: 
:Created: None
:Updated: None



Outputs
-----------

========  ========  ===============
name      type      documentation
========  ========  ===============
bwa       FastaBwa
samtools  FastaFai
dict      FastDict
========  ========  ===============


Embedded Tools
***************

===============================  =========================================
BWA-Index                        ``bwaIndex/v0.7.15``
SamTools: faidx                  ``SamToolsIndex/1.9.0``
GATK4: CreateSequenceDictionary  ``Gatk4CreateSequenceDictionary/4.1.3.0``
===============================  =========================================



Additional configuration (inputs)
---------------------------------

=========  ======  ===============
name       type    documentation
=========  ======  ===============
reference  Fasta
=========  ======  ===============


