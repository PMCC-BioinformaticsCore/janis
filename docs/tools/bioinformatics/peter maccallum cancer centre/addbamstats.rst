:orphan:

Add Bam Statistics to Vcf
=======================================

``addBamStats`` · *1 contributor · 1 version*

usage: add_bam_stats.py [-h] -i I -o O --type {germline,somatic}
                        [--mpileup MPILEUP] [--normal_mpileup NORMAL_MPILEUP]
                        [--tumor_mpileup TUMOR_MPILEUP]
                        [--normal_id NORMAL_ID] [--tumor_id TUMOR_ID]

Get stats from bam file and write to vcf

required arguments:
  -i I                  input vcf
  -o O                  output vcf
  --type {germline,somatic}
                        must be either germline or somatic
  --mpileup MPILEUP     mpileup file extracted from bam file
  --normal_mpileup NORMAL_MPILEUP
                        mpileup file extracted from the normal sample bam,
                        required if input is somatic vcf
  --tumor_mpileup TUMOR_MPILEUP
                        mpileup file extracted from the tumor sample, required
                        if input is somatic vcf
  --normal_id NORMAL_ID
                        Normal sample id, required if input is somatic vcf
  --tumor_id TUMOR_ID   Tumor sample id, required if input is somatic vcf

optional arguments:
  -h, --help            show this help message and exit
        


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.addbamstats.versions import AddBamStats_0_0_7

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "addbamstats_step",
           AddBamStats_0_0_7(
               inputVcf=None,
               type=None,
           )
       )
       wf.output("out", source=addbamstats_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for addBamStats:

.. code-block:: bash

   # user inputs
   janis inputs addBamStats > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inputVcf: inputVcf.vcf
       type: <value>




5. Run addBamStats with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       addBamStats





Information
------------

:ID: ``addBamStats``
:URL: `https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/vcf_utils <https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/vcf_utils>`_
:Versions: 0.0.7
:Container: michaelfranklin/pmacutil:0.0.7
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-05-20 00:00:00
:Updated: 2020-05-20 00:00:00


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============


Additional configuration (inputs)
---------------------------------

==============  ==================  ================  ==========  ===================================================================================
name            type                prefix            position    documentation
==============  ==================  ================  ==========  ===================================================================================
inputVcf        VCF                 -i                            input vcf
type            String              --type                        must be either germline or somatic
mpileup         Optional<File>      --mpileup                     mpileup file extracted from bam file
normalMpileup   Optional<File>      --normal_mpileup              mpileup file extracted from the normal sample bam, required if input is somatic vcf
tumorMpileup    Optional<File>      --tumor_mpileup               mpileup file extracted from the tumor sample bam, required if input is somatic vcf
normalID        Optional<String>    --normal_id                   normal sample id, required if input is somatic vcf
tumorID         Optional<String>    --tumor_id                    tumor sample id, required if input is somatic vcf
outputFilename  Optional<Filename>  -o                            output vcf name
==============  ==================  ================  ==========  ===================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task addBamStats {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File? mpileup
       File? normalMpileup
       File? tumorMpileup
       String? normalID
       String? tumorID
       File inputVcf
       String? outputFilename
       String type
     }
     command <<<
       set -e
       add_bam_stats.py \
         ~{if defined(mpileup) then ("--mpileup '" + mpileup + "'") else ""} \
         ~{if defined(normalMpileup) then ("--normal_mpileup '" + normalMpileup + "'") else ""} \
         ~{if defined(tumorMpileup) then ("--tumor_mpileup '" + tumorMpileup + "'") else ""} \
         ~{if defined(normalID) then ("--normal_id '" + normalID + "'") else ""} \
         ~{if defined(tumorID) then ("--tumor_id '" + tumorID + "'") else ""} \
         -i '~{inputVcf}' \
         -o '~{select_first([outputFilename, "generated.addbamstats.vcf"])}' \
         --type '~{type}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "michaelfranklin/pmacutil:0.0.7"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.addbamstats.vcf"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Add Bam Statistics to Vcf
   doc: |-
     usage: add_bam_stats.py [-h] -i I -o O --type {germline,somatic}
                             [--mpileup MPILEUP] [--normal_mpileup NORMAL_MPILEUP]
                             [--tumor_mpileup TUMOR_MPILEUP]
                             [--normal_id NORMAL_ID] [--tumor_id TUMOR_ID]

     Get stats from bam file and write to vcf

     required arguments:
       -i I                  input vcf
       -o O                  output vcf
       --type {germline,somatic}
                             must be either germline or somatic
       --mpileup MPILEUP     mpileup file extracted from bam file
       --normal_mpileup NORMAL_MPILEUP
                             mpileup file extracted from the normal sample bam,
                             required if input is somatic vcf
       --tumor_mpileup TUMOR_MPILEUP
                             mpileup file extracted from the tumor sample, required
                             if input is somatic vcf
       --normal_id NORMAL_ID
                             Normal sample id, required if input is somatic vcf
       --tumor_id TUMOR_ID   Tumor sample id, required if input is somatic vcf

     optional arguments:
       -h, --help            show this help message and exit
          

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: michaelfranklin/pmacutil:0.0.7

   inputs:
   - id: mpileup
     label: mpileup
     doc: mpileup file extracted from bam file
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --mpileup
   - id: normalMpileup
     label: normalMpileup
     doc: |-
       mpileup file extracted from the normal sample bam, required if input is somatic vcf
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --normal_mpileup
   - id: tumorMpileup
     label: tumorMpileup
     doc: |-
       mpileup file extracted from the tumor sample bam, required if input is somatic vcf
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --tumor_mpileup
   - id: normalID
     label: normalID
     doc: normal sample id, required if input is somatic vcf
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --normal_id
   - id: tumorID
     label: tumorID
     doc: tumor sample id, required if input is somatic vcf
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --tumor_id
   - id: inputVcf
     label: inputVcf
     doc: input vcf
     type: File
     inputBinding:
       prefix: -i
   - id: outputFilename
     label: outputFilename
     doc: output vcf name
     type:
     - string
     - 'null'
     default: generated.addbamstats.vcf
     inputBinding:
       prefix: -o
   - id: type
     label: type
     doc: must be either germline or somatic
     type: string
     inputBinding:
       prefix: --type

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.addbamstats.vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: add_bam_stats.py
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: addBamStats


