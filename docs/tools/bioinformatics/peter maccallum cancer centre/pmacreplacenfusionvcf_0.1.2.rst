:orphan:

Peter Mac: replace_N_fusion_vcf.py
==========================================================

``PMacReplaceNFusionVcf`` · *1 contributor · 1 version*


Copy the REF base over to ALT colum for breakends sv
Usage: replace_N_fusion_vcf.py --input $input_vcf > $output_vcf



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.replacenfusionvcf.versions import ReplaceNFusionVcf_0_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "pmacreplacenfusionvcf_step",
           ReplaceNFusionVcf_0_1_2(
               vcf=None,
           )
       )
       wf.output("out", source=pmacreplacenfusionvcf_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for PMacReplaceNFusionVcf:

.. code-block:: bash

   # user inputs
   janis inputs PMacReplaceNFusionVcf > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf




5. Run PMacReplaceNFusionVcf with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       PMacReplaceNFusionVcf

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          PMacReplaceNFusionVcf

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``PMacReplaceNFusionVcf``
:URL: *No URL to the documentation was provided*
:Versions: 0.1.2
:Container: rlupat/pmacutil:latest
:Authors: Jiaan Yu
:Citations: None
:Created: 2022-01-05
:Updated: 2022-01-05


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============


Additional configuration (inputs)
---------------------------------

==============  ==================  ========  ==========  ===============
name            type                prefix      position  documentation
==============  ==================  ========  ==========  ===============
vcf             VCF                 --input            2
outputFilename  Optional<Filename>                     5
==============  ==================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task PMacReplaceNFusionVcf {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       File vcf
       String? outputFilename
     }

     command <<<
       set -e
       replace_N_fusion_vcf.py \
         --input '~{vcf}' \
         > \
         '~{select_first([outputFilename, "~{basename(vcf, ".vcf")}.vcf"])}'
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "rlupat/pmacutil:latest"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }

     output {
       File out = select_first([outputFilename, "~{basename(vcf, ".vcf")}.vcf"])
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'Peter Mac: replace_N_fusion_vcf.py'

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: rlupat/pmacutil:latest

   inputs:
   - id: vcf
     label: vcf
     type: File
     inputBinding:
       prefix: --input
       position: 2
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.vcf
     inputBinding:
       position: 5
       valueFrom: $(inputs.vcf.basename.replace(/.vcf$/, "")).vcf

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: $(inputs.vcf.basename.replace(/.vcf$/, "")).vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - replace_N_fusion_vcf.py
   arguments:
   - position: 4
     valueFrom: '>'
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: PMacReplaceNFusionVcf


