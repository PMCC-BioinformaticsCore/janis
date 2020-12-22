:orphan:

Awk
=========

``awk`` · *0 contributors · 1 version*

run an awk script


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.awk import Awk

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "awk_step",
           Awk(
               script=None,
               input_files=None,
           )
       )
       wf.output("out", source=awk_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for awk:

.. code-block:: bash

   # user inputs
   janis inputs awk > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       input_files:
       - input_files_0
       - input_files_1
       script: script




5. Run awk with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       awk





Information
------------

:ID: ``awk``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956
:Authors: 
:Citations: None
:Created: None
:Updated: None


Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     stdout<File>
======  ============  ===============


Additional configuration (inputs)
---------------------------------

===========  ===========  ========  ==========  ===============
name         type         prefix      position  documentation
===========  ===========  ========  ==========  ===============
script       File         -f                 1
input_files  Array<File>                     2
===========  ===========  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task awk {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File script
       Array[File] input_files
     }
     command <<<
       set -e
       awk \
         -f '~{script}' \
         ~{if length(input_files) > 0 then "'" + sep("' '", input_files) + "'" else ""}
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
       File out = stdout()
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Awk
   doc: run an awk script

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

   inputs:
   - id: script
     label: script
     type: File
     inputBinding:
       prefix: -f
       position: 1
   - id: input_files
     label: input_files
     type:
       type: array
       items: File
     inputBinding:
       position: 2

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: awk
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: awk


