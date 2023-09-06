:orphan:

GatherFilesForMultiqc
=====================

``GatherFilesForMultiqc`` · *1 contributor · 1 version*

Gather Files for MultiQC
        


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.gatherfilesformultiqc import GatherFilesForMultiqc

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatherfilesformultiqc_step",
           GatherFilesForMultiqc(
               inp_files=None,
               inp_files2=None,
           )
       )
       wf.output("out", source=gatherfilesformultiqc_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for GatherFilesForMultiqc:

.. code-block:: bash

   # user inputs
   janis inputs GatherFilesForMultiqc > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inp_files:
       - inp_files_0
       - inp_files_1
       inp_files2:
       - inp_files2_0
       - inp_files2_1




5. Run GatherFilesForMultiqc with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       GatherFilesForMultiqc

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          GatherFilesForMultiqc

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``GatherFilesForMultiqc``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Container: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956
:Authors: Jiaan Yu
:Citations: None
:Created: 2021-11-01 00:00:00
:Updated: 2021-11-19 00:00:00


Outputs
-----------

======  =========  ===============
name    type       documentation
======  =========  ===============
out     Directory
======  =========  ===============


Additional configuration (inputs)
---------------------------------

==========  ================  ========  ==========  ===============
name        type              prefix      position  documentation
==========  ================  ========  ==========  ===============
inp_files   Array<File>                          4
inp_files2  Array<File>                          5
output_dir  Optional<String>                     8
==========  ================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task GatherFilesForMultiqc {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       Array[File] inp_files
       Array[File] inp_files2
       String? output_dir
     }

     command <<<
       set -e
        \
         mkdir \
         ~{select_first([output_dir, "output_dir"])} \
         ; \
         cp \
         ~{if length(inp_files) > 0 then "'" + sep("' '", inp_files) + "'" else ""} \
         ~{if length(inp_files2) > 0 then "'" + sep("' '", inp_files2) + "'" else ""} \
         ~{if defined(select_first([output_dir, "output_dir"])) then ("'" + select_first([output_dir, "output_dir"]) + "'") else ""}
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }

     output {
       Directory out = select_first([output_dir, "output_dir"])
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: GatherFilesForMultiqc

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

   inputs:
   - id: inp_files
     label: inp_files
     type:
       type: array
       items: File
     inputBinding:
       position: 4
   - id: inp_files2
     label: inp_files2
     type:
       type: array
       items: File
     inputBinding:
       position: 5
   - id: output_dir
     label: output_dir
     type: string
     default: output_dir
     inputBinding:
       position: 8

   outputs:
   - id: out
     label: out
     type: Directory
     outputBinding:
       glob: |-
         $((inputs.output_dir ? inputs.output_dir : "generated" != null) ? inputs.output_dir ? inputs.output_dir : "generated" : "output_dir")
       loadContents: false
   stdout: _stdout
   stderr: _stderr
   arguments:
   - position: 0
     valueFrom: mkdir
     shellQuote: false
   - position: 1
     valueFrom: $(inputs.output_dir)
     shellQuote: false
   - position: 2
     valueFrom: ;
     shellQuote: false
   - position: 3
     valueFrom: cp
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: GatherFilesForMultiqc


