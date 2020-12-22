:orphan:

Tar (archive)
===================

``Tar`` · *0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-unix>`_


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.tar import Tar

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "tar_step",
           Tar(
               files=None,
           )
       )
       wf.output("out", source=tar_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Tar:

.. code-block:: bash

   # user inputs
   janis inputs Tar > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       files:
       - files_0
       - files_1




5. Run Tar with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Tar





Information
------------

:ID: ``Tar``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956
:Authors: 
:Citations: None
:Created: None
:Updated: None


Outputs
-----------

======  =======  ===============
name    type     documentation
======  =======  ===============
out     TarFile
======  =======  ===============


Additional configuration (inputs)
---------------------------------

==============  =====================  ========  ==========  ===============
name            type                   prefix      position  documentation
==============  =====================  ========  ==========  ===============
files           Array<File>                               2
files2          Optional<Array<File>>                     3
outputFilename  Optional<Filename>                        1
==============  =====================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Tar {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[File] files
       Array[File]? files2
       String? outputFilename
     }
     command <<<
       set -e
       tar cvf \
         '~{select_first([outputFilename, "generated.tar"])}' \
         ~{if length(files) > 0 then "'" + sep("' '", files) + "'" else ""} \
         ~{if (defined(files2) && length(select_first([files2])) > 0) then "'" + sep("' '", select_first([files2])) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.tar"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Tar (archive)

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

   inputs:
   - id: files
     label: files
     type:
       type: array
       items: File
     inputBinding:
       position: 2
   - id: files2
     label: files2
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       position: 3
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.tar
     inputBinding:
       position: 1

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.tar
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - tar
   - cvf
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Tar


