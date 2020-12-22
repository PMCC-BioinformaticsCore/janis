:orphan:

Extract Strelka Somatic AD DP
=========================================================

``extractStrelkaSomaticADDP`` · *1 contributor · 2 versions*


 - Extract and calculate AD and AF value for each variant (both SNVs and INDELs)
Based on https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#somatic
        


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.extractstrelkasomaticaddp.versions import ExtractStrelkaSomaticADDP_0_1_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "extractstrelkasomaticaddp_step",
           ExtractStrelkaSomaticADDP_0_1_1(
               vcf=None,
           )
       )
       wf.output("out", source=extractstrelkasomaticaddp_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for extractStrelkaSomaticADDP:

.. code-block:: bash

   # user inputs
   janis inputs extractStrelkaSomaticADDP > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf




5. Run extractStrelkaSomaticADDP with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       extractStrelkaSomaticADDP





Information
------------

:ID: ``extractStrelkaSomaticADDP``
:URL: `https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/vcf_utils <https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/vcf_utils>`_
:Versions: 0.1.1, 0.1.0
:Container: michaelfranklin/pmacutil:0.1.1
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-07-27 00:00:00
:Updated: 2020-07-27 00:00:00


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
name            type                prefix    position    documentation
==============  ==================  ========  ==========  ===============
vcf             VCF                 -i                    input vcf
outputFilename  Optional<Filename>  -o                    output vcf
==============  ==================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task extractStrelkaSomaticADDP {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
       String? outputFilename
     }
     command <<<
       set -e
       extract_strelka_somatic_DP_AF.py \
         -i '~{vcf}' \
         -o '~{select_first([outputFilename, "generated.vcf"])}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "michaelfranklin/pmacutil:0.1.1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.vcf"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Extract Strelka Somatic AD DP
   doc: |2-

      - Extract and calculate AD and AF value for each variant (both SNVs and INDELs)
     Based on https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#somatic
          

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: michaelfranklin/pmacutil:0.1.1

   inputs:
   - id: vcf
     label: vcf
     doc: input vcf
     type: File
     inputBinding:
       prefix: -i
   - id: outputFilename
     label: outputFilename
     doc: output vcf
     type:
     - string
     - 'null'
     default: generated.vcf
     inputBinding:
       prefix: -o

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: extract_strelka_somatic_DP_AF.py
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: extractStrelkaSomaticADDP


