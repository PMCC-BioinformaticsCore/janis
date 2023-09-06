:orphan:

Peter Mac: MegaFusion.py
=========================================

``PMacMegaFusion`` · *1 contributor · 1 version*


Convert RNA fusion files to SV VCF. MegaFusion accepts a fusion transcript 
file produced from any tool (such as Arriba), as well as a JSON file, 
specifying which columns to put in the output vcf file.
Usage: MegaFusion.py --sample $sample --json $json         --fusion $fusion --contig $contig > $output_vcf



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.megafusion.versions import MegaFusion_0_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "pmacmegafusion_step",
           MegaFusion_0_1_2(
               json=None,
               fusion=None,
               toolVersion=None,
           )
       )
       wf.output("out", source=pmacmegafusion_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for PMacMegaFusion:

.. code-block:: bash

   # user inputs
   janis inputs PMacMegaFusion > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fusion: fusion.tsv
       json: json
       toolVersion: <value>




5. Run PMacMegaFusion with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       PMacMegaFusion

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          PMacMegaFusion

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``PMacMegaFusion``
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

==============  ==================  ==============  ==========  ===============
name            type                prefix            position  documentation
==============  ==================  ==============  ==========  ===============
json            File                --json                   2
fusion          tsv                 --fusion                 2
toolVersion     String              --tool_version           2
sample          Optional<String>    --sample                 2
contig          Optional<File>      --contig                 2
outputFilename  Optional<Filename>                           5
==============  ==================  ==============  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task PMacMegaFusion {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       String? sample
       File json
       File fusion
       String toolVersion
       File? contig
       String? outputFilename
     }

     command <<<
       set -e
       MegaFusion.py \
         ~{if defined(sample) then ("--sample '" + sample + "'") else ""} \
         --json '~{json}' \
         --fusion '~{fusion}' \
         --tool_version '~{toolVersion}' \
         ~{if defined(contig) then ("--contig '" + contig + "'") else ""} \
         > \
         '~{select_first([outputFilename, "~{basename(fusion, ".tsv")}.vcf"])}'
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
       File out = select_first([outputFilename, "~{basename(fusion, ".tsv")}.vcf"])
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'Peter Mac: MegaFusion.py'

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: rlupat/pmacutil:latest

   inputs:
   - id: sample
     label: sample
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sample
       position: 2
   - id: json
     label: json
     type: File
     inputBinding:
       prefix: --json
       position: 2
   - id: fusion
     label: fusion
     type: File
     inputBinding:
       prefix: --fusion
       position: 2
   - id: toolVersion
     label: toolVersion
     type: string
     inputBinding:
       prefix: --tool_version
       position: 2
   - id: contig
     label: contig
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --contig
       position: 2
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.vcf
     inputBinding:
       position: 5
       valueFrom: $(inputs.fusion.basename.replace(/.tsv$/, "")).vcf

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: $(inputs.fusion.basename.replace(/.tsv$/, "")).vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - MegaFusion.py
   arguments:
   - position: 4
     valueFrom: '>'
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: PMacMegaFusion


