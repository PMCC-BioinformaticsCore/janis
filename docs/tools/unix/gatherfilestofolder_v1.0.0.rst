:orphan:

GatherFilesToFolder
===================

``GatherFilesToFolder`` · *0 contributors · 1 version*

gather files to a folder


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.gatherfilestofolder import GatherFilesToFolder

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatherfilestofolder_step",
           GatherFilesToFolder(
               inp_files=None,
           )
       )
       wf.output("out", source=gatherfilestofolder_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for GatherFilesToFolder:

.. code-block:: bash

   # user inputs
   janis inputs GatherFilesToFolder > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inp_files:
       - inp_files_0
       - inp_files_1




5. Run GatherFilesToFolder with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       GatherFilesToFolder

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          GatherFilesToFolder

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``GatherFilesToFolder``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956
:Authors: 
:Citations: None
:Created: None
:Updated: None


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
output_dir  Optional<String>                     8
==========  ================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task GatherFilesToFolder {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       Array[File] inp_files
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
   label: GatherFilesToFolder

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
   id: GatherFilesToFolder


