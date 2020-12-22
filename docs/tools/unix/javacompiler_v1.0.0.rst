:orphan:

Java compiler
============================

``javacompiler`` · *0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-unix>`_


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.compile import Compile

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "javacompiler_step",
           Compile(
               file=None,
           )
       )
       wf.output("out", source=javacompiler_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for javacompiler:

.. code-block:: bash

   # user inputs
   janis inputs javacompiler > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       file: file




5. Run javacompiler with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       javacompiler





Information
------------

:ID: ``javacompiler``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: openjdk:8
:Authors: 
:Citations: None
:Created: None
:Updated: None


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
======  ======  ===============


Additional configuration (inputs)
---------------------------------

======  ======  ========  ==========  ===============
name    type    prefix      position  documentation
======  ======  ========  ==========  ===============
file    File                       1
======  ======  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task javacompiler {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File file
     }
     command <<<
       set -e
       javac \
         -d '.' \
         '~{file}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "openjdk:8"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = glob("*.class")[0]
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Java compiler

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: openjdk:8

   inputs:
   - id: file
     label: file
     type: File
     inputBinding:
       position: 1

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: '*.class'
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: javac
   arguments:
   - prefix: -d
     position: 0
     valueFrom: .

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: javacompiler


