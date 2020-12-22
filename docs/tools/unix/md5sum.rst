:orphan:

MD5 Sum
================

``md5sum`` · *0 contributors · 1 version*

Compute the MD5 message digest of the given file.


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.md5sum import MD5Sum

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "md5sum_step",
           MD5Sum(
               input_file=None,
           )
       )
       wf.output("out", source=md5sum_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for md5sum:

.. code-block:: bash

   # user inputs
   janis inputs md5sum > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       input_file: input_file




5. Run md5sum with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       md5sum





Information
------------

:ID: ``md5sum``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956
:Authors: 
:Citations: None
:Created: None
:Updated: 2020-06-09 00:00:00


Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     stdout<File>
======  ============  ===============


Additional configuration (inputs)
---------------------------------

==========  ======  ========  ==========  ===============
name        type    prefix      position  documentation
==========  ======  ========  ==========  ===============
input_file  File                       1
==========  ======  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task md5sum {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File input_file
     }
     command <<<
       set -e
       md5sum \
         '~{input_file}' \
         | awk '{print $1}'
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
   label: MD5 Sum
   doc: Compute the MD5 message digest of the given file.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

   inputs:
   - id: input_file
     label: input_file
     type: File
     inputBinding:
       position: 1

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: md5sum
   arguments:
   - position: 2
     valueFrom: "| awk '{print $1}'"
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: md5sum


