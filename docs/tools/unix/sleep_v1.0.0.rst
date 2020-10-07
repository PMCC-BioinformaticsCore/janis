:orphan:

Sleep
=============

``sleep`` · *0 contributors · 1 version*

sleep for the given number of seconds


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.sleep import Sleep

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "sleep_step",
           Sleep(
               time=None,
           )
       )
       wf.output("out", source=sleep_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for sleep:

.. code-block:: bash

   # user inputs
   janis inputs sleep > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       time: 0




5. Run sleep with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       sleep





Information
------------

:ID: ``sleep``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu:latest
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

======  =======  ========  ==========  ===============
name    type     prefix      position  documentation
======  =======  ========  ==========  ===============
time    Integer                     1
======  =======  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task sleep {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Int time
     }
     command <<<
       set -e
       sleep \
         ~{time}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "ubuntu:latest"
       duration: select_first([runtime_seconds, (time + 30), 86400])
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
   cwlVersion: v1.0
   label: Sleep
   doc: sleep for the given number of seconds

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ubuntu:latest

   inputs:
   - id: time
     label: time
     type: int
     inputBinding:
       position: 1

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: sleep
   arguments: []
   id: sleep


