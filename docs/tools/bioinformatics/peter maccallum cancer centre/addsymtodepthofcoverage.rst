:orphan:

Add Sym to DepthOfCoverage
====================================================

``addSymToDepthOfCoverage`` · *1 contributor · 1 version*

usage: add_sym_to_DepthOfCoverage.py [-h] -i INPUT -o OUTPUT -bed BED

Performance summary of bam

optional arguments:
  -h, --help  show this help message and exit
  -i INPUT    Gatk3 DepthOfCoverage interval_summary output
  -o OUTPUT   Output file name
  -bed BED    Annotated bed file
        


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.addsymtodepthofcoverage.versions import AddSymToDepthOfCoverage_0_0_7

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "addsymtodepthofcoverage_step",
           AddSymToDepthOfCoverage_0_0_7(
               inputFile=None,
               bed=None,
           )
       )
       wf.output("out", source=addsymtodepthofcoverage_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for addSymToDepthOfCoverage:

.. code-block:: bash

   # user inputs
   janis inputs addSymToDepthOfCoverage > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bed: bed.bed
       inputFile: inputFile




5. Run addSymToDepthOfCoverage with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       addSymToDepthOfCoverage





Information
------------

:ID: ``addSymToDepthOfCoverage``
:URL: `https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/performance <https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/performance>`_
:Versions: 0.0.7
:Container: michaelfranklin/pmacutil:0.0.7
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-04-09 00:00:00
:Updated: 2020-04-09 00:00:00


Outputs
-----------

======  ========  ===============
name    type      documentation
======  ========  ===============
out     TextFile
======  ========  ===============


Additional configuration (inputs)
---------------------------------

==============  ==================  ========  ==========  =============================================
name            type                prefix    position    documentation
==============  ==================  ========  ==========  =============================================
inputFile       File                -i                    Gatk3 DepthOfCoverage interval_summary output
bed             bed                 -bed                  Annotated bed file
outputFilename  Optional<Filename>  -o                    Output file name
==============  ==================  ========  ==========  =============================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task addSymToDepthOfCoverage {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File inputFile
       String? outputFilename
       File bed
     }
     command <<<
       set -e
       add_sym_to_DepthOfCoverage.py \
         -i '~{inputFile}' \
         -o '~{select_first([outputFilename, "generated.txt"])}' \
         -bed '~{bed}'
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
       File out = select_first([outputFilename, "generated.txt"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Add Sym to DepthOfCoverage
   doc: |-
     usage: add_sym_to_DepthOfCoverage.py [-h] -i INPUT -o OUTPUT -bed BED

     Performance summary of bam

     optional arguments:
       -h, --help  show this help message and exit
       -i INPUT    Gatk3 DepthOfCoverage interval_summary output
       -o OUTPUT   Output file name
       -bed BED    Annotated bed file
          

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: michaelfranklin/pmacutil:0.0.7

   inputs:
   - id: inputFile
     label: inputFile
     doc: Gatk3 DepthOfCoverage interval_summary output
     type: File
     inputBinding:
       prefix: -i
   - id: outputFilename
     label: outputFilename
     doc: Output file name
     type:
     - string
     - 'null'
     default: generated.txt
     inputBinding:
       prefix: -o
   - id: bed
     label: bed
     doc: Annotated bed file
     type: File
     inputBinding:
       prefix: -bed

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.txt
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: add_sym_to_DepthOfCoverage.py
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: addSymToDepthOfCoverage


