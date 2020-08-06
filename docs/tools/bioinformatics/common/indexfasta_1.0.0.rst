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
SamTools: faidx                  ``SamToolsFaidx/1.9.0``
GATK4: CreateSequenceDictionary  ``Gatk4CreateSequenceDictionary/4.1.3.0``
===============================  =========================================



Additional configuration (inputs)
---------------------------------

=========  ======  ===============
name       type    documentation
=========  ======  ===============
reference  Fasta
=========  ======  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   import "tools/bwaIndex_v0_7_15.wdl" as B
   import "tools/SamToolsFaidx_1_9_0.wdl" as S
   import "tools/Gatk4CreateSequenceDictionary_4_1_3_0.wdl" as G

   workflow IndexFasta {
     input {
       File reference
     }
     call B.bwaIndex as create_bwa {
       input:
         reference=reference
     }
     call S.SamToolsFaidx as create_samtools {
       input:
         reference=reference
     }
     call G.Gatk4CreateSequenceDictionary as create_dict {
       input:
         reference=reference
     }
     output {
       File bwa = create_bwa.out
       File bwa_amb = create_bwa.out_amb
       File bwa_ann = create_bwa.out_ann
       File bwa_bwt = create_bwa.out_bwt
       File bwa_pac = create_bwa.out_pac
       File bwa_sa = create_bwa.out_sa
       File samtools = create_samtools.out
       File samtools_fai = create_samtools.out_fai
       File dict = create_dict.out
       File dict_dict = create_dict.out_dict
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: Workflow
   cwlVersion: v1.0
   label: Index Fasta reference

   requirements:
   - class: InlineJavascriptRequirement
   - class: StepInputExpressionRequirement

   inputs:
   - id: reference
     type: File

   outputs:
   - id: bwa
     type: File
     secondaryFiles:
     - .amb
     - .ann
     - .bwt
     - .pac
     - .sa
     outputSource: create_bwa/out
   - id: samtools
     type: File
     secondaryFiles:
     - .fai
     outputSource: create_samtools/out
   - id: dict
     type: File
     secondaryFiles:
     - ^.dict
     outputSource: create_dict/out

   steps:
   - id: create_bwa
     label: BWA-Index
     in:
     - id: reference
       source: reference
     run: tools/bwaIndex_v0_7_15.cwl
     out:
     - id: out
   - id: create_samtools
     label: 'SamTools: faidx'
     in:
     - id: reference
       source: reference
     run: tools/SamToolsFaidx_1_9_0.cwl
     out:
     - id: out
   - id: create_dict
     label: 'GATK4: CreateSequenceDictionary'
     in:
     - id: reference
       source: reference
     run: tools/Gatk4CreateSequenceDictionary_4_1_3_0.cwl
     out:
     - id: out
   id: IndexFasta

