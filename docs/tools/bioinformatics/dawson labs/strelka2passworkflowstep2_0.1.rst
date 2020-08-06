:orphan:

Strelka 2Pass analysis step 2
=========================================================

*1 contributor · 1 version*

This is the second step for joint somatic variant calling
        based on a 2pass analysis common in RNASeq.

        It runs strelka2 again with the variants found in all of the other samples as input to be forced to genotype these.

        It also normalises and indexes the output vcfs


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.dawson.workflows.strelka2passanalysisstep2 import Strelka2PassWorkflowStep2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "strelka2passworkflowstep2_step",
           Strelka2PassWorkflowStep2(
               normalBam=None,
               tumorBam=None,
               reference=None,
               indelCandidates=None,
               strelkaSNVs=None,
           )
       )
       wf.output("indels", source=strelka2passworkflowstep2_step.indels)
       wf.output("snvs", source=strelka2passworkflowstep2_step.snvs)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Strelka2PassWorkflowStep2:

.. code-block:: bash

   # user inputs
   janis inputs Strelka2PassWorkflowStep2 > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       indelCandidates:
       - indelCandidates_0.vcf.gz
       - indelCandidates_1.vcf.gz
       normalBam: normalBam.cram
       reference: reference.fasta
       strelkaSNVs:
       - strelkaSNVs_0.vcf.gz
       - strelkaSNVs_1.vcf.gz
       tumorBam: tumorBam.cram




5. Run Strelka2PassWorkflowStep2 with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Strelka2PassWorkflowStep2





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``Strelka2PassWorkflowStep2``
:URL: *No URL to the documentation was provided*
:Versions: 0.1
:Authors: Sebastian Hollizeck
:Citations: 
:Created: 2019-10-11
:Updated: 2020-08-04



Outputs
-----------

======  ====================  ===============
name    type                  documentation
======  ====================  ===============
indels  CompressedIndexedVCF
snvs    CompressedIndexedVCF
======  ====================  ===============


Embedded Tools
***************

===================  ===============================
Strelka (Somatic)    ``strelka_somatic_cram/2.9.10``
BCFTools: Normalize  ``bcftoolsNorm/v1.9``
BCFTools: Index      ``bcftoolsIndex/v1.9``
===================  ===============================



Additional configuration (inputs)
---------------------------------

===============  ===========================  ===============
name             type                         documentation
===============  ===========================  ===============
normalBam        CramPair
tumorBam         CramPair
reference        FastaWithIndexes
indelCandidates  Array<CompressedIndexedVCF>
strelkaSNVs      Array<CompressedIndexedVCF>
callRegions      Optional<BedTABIX>
exome            Optional<Boolean>
configStrelka    Optional<File>
===============  ===========================  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   import "tools/strelka_somatic_cram_2_9_10.wdl" as S
   import "tools/bcftoolsNorm_v1_9.wdl" as B
   import "tools/bcftoolsIndex_v1_9.wdl" as B2

   workflow Strelka2PassWorkflowStep2 {
     input {
       File normalBam
       File normalBam_crai
       File tumorBam
       File tumorBam_crai
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       File? callRegions
       File? callRegions_tbi
       Boolean? exome = false
       File? configStrelka
       Array[File] indelCandidates
       Array[File] indelCandidates_tbi
       Array[File] strelkaSNVs
       Array[File] strelkaSNVs_tbi
     }
     call S.strelka_somatic_cram as strelka2pass {
       input:
         normalBam=normalBam,
         normalBam_crai=normalBam_crai,
         tumorBam=tumorBam,
         tumorBam_crai=tumorBam_crai,
         reference=reference,
         reference_fai=reference_fai,
         reference_amb=reference_amb,
         reference_ann=reference_ann,
         reference_bwt=reference_bwt,
         reference_pac=reference_pac,
         reference_sa=reference_sa,
         reference_dict=reference_dict,
         config=configStrelka,
         indelCandidates=indelCandidates,
         indelCandidates_tbi=indelCandidates_tbi,
         forcedgt=strelkaSNVs,
         forcedgt_tbi=strelkaSNVs_tbi,
         exome=select_first([exome, false]),
         callRegions=callRegions,
         callRegions_tbi=callRegions_tbi
     }
     call B.bcftoolsNorm as normaliseSNVs {
       input:
         vcf=strelka2pass.snvs,
         reference=reference,
         reference_fai=reference_fai,
         reference_amb=reference_amb,
         reference_ann=reference_ann,
         reference_bwt=reference_bwt,
         reference_pac=reference_pac,
         reference_sa=reference_sa,
         reference_dict=reference_dict
     }
     call B2.bcftoolsIndex as indexSNVs {
       input:
         vcf=normaliseSNVs.out
     }
     call B.bcftoolsNorm as normaliseINDELs {
       input:
         vcf=strelka2pass.indels,
         reference=reference,
         reference_fai=reference_fai,
         reference_amb=reference_amb,
         reference_ann=reference_ann,
         reference_bwt=reference_bwt,
         reference_pac=reference_pac,
         reference_sa=reference_sa,
         reference_dict=reference_dict
     }
     call B2.bcftoolsIndex as indexINDELs {
       input:
         vcf=normaliseINDELs.out
     }
     output {
       File indels = indexINDELs.out
       File indels_tbi = indexINDELs.out_tbi
       File snvs = indexSNVs.out
       File snvs_tbi = indexSNVs.out_tbi
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: Workflow
   cwlVersion: v1.0
   label: Strelka 2Pass analysis step 2
   doc: |-
     This is the second step for joint somatic variant calling
             based on a 2pass analysis common in RNASeq.

             It runs strelka2 again with the variants found in all of the other samples as input to be forced to genotype these.

             It also normalises and indexes the output vcfs

   requirements:
   - class: InlineJavascriptRequirement
   - class: StepInputExpressionRequirement

   inputs:
   - id: normalBam
     type: File
     secondaryFiles:
     - .crai
   - id: tumorBam
     type: File
     secondaryFiles:
     - .crai
   - id: reference
     type: File
     secondaryFiles:
     - .fai
     - .amb
     - .ann
     - .bwt
     - .pac
     - .sa
     - ^.dict
   - id: callRegions
     type:
     - File
     - 'null'
     secondaryFiles:
     - .tbi
   - id: exome
     type: boolean
     default: false
   - id: configStrelka
     type:
     - File
     - 'null'
   - id: indelCandidates
     type:
       type: array
       items: File
     secondaryFiles:
     - .tbi
   - id: strelkaSNVs
     type:
       type: array
       items: File
     secondaryFiles:
     - .tbi

   outputs:
   - id: indels
     type: File
     secondaryFiles:
     - .tbi
     outputSource: indexINDELs/out
   - id: snvs
     type: File
     secondaryFiles:
     - .tbi
     outputSource: indexSNVs/out

   steps:
   - id: strelka2pass
     label: Strelka (Somatic)
     in:
     - id: normalBam
       source: normalBam
     - id: tumorBam
       source: tumorBam
     - id: reference
       source: reference
     - id: config
       source: configStrelka
     - id: indelCandidates
       source: indelCandidates
     - id: forcedgt
       source: strelkaSNVs
     - id: exome
       source: exome
     - id: callRegions
       source: callRegions
     run: tools/strelka_somatic_cram_2_9_10.cwl
     out:
     - id: configPickle
     - id: script
     - id: stats
     - id: indels
     - id: snvs
   - id: normaliseSNVs
     label: 'BCFTools: Normalize'
     in:
     - id: vcf
       source: strelka2pass/snvs
     - id: reference
       source: reference
     run: tools/bcftoolsNorm_v1_9.cwl
     out:
     - id: out
   - id: indexSNVs
     label: 'BCFTools: Index'
     in:
     - id: vcf
       source: normaliseSNVs/out
     run: tools/bcftoolsIndex_v1_9.cwl
     out:
     - id: out
   - id: normaliseINDELs
     label: 'BCFTools: Normalize'
     in:
     - id: vcf
       source: strelka2pass/indels
     - id: reference
       source: reference
     run: tools/bcftoolsNorm_v1_9.cwl
     out:
     - id: out
   - id: indexINDELs
     label: 'BCFTools: Index'
     in:
     - id: vcf
       source: normaliseINDELs/out
     run: tools/bcftoolsIndex_v1_9.cwl
     out:
     - id: out
   id: Strelka2PassWorkflowStep2

