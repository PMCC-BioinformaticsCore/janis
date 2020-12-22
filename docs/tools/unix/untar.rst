:orphan:

Tar (unarchive)
=======================

``untar`` · *0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-unix>`_


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.untar import Untar

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "untar_step",
           Untar(
               tarfile=None,
           )
       )
       wf.output("out", source=untar_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for untar:

.. code-block:: bash

   # user inputs
   janis inputs untar > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       tarfile: tarfile.tar




5. Run untar with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       untar





Information
------------

:ID: ``untar``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956
:Authors: 
:Citations: None
:Created: None
:Updated: None


Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     Array<File>
======  ===========  ===============


Additional configuration (inputs)
---------------------------------

=======  =======  ========  ==========  ===============
name     type     prefix      position  documentation
=======  =======  ========  ==========  ===============
tarfile  TarFile                     0
=======  =======  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task untar {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File tarfile
     }
     command <<<
       set -e
       tar xf \
         '~{tarfile}'
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
       Array[File] out = glob("*.java")
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Tar (unarchive)

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

   inputs:
   - id: tarfile
     label: tarfile
     type: File
     inputBinding:
       position: 0

   outputs:
   - id: out
     label: out
     type:
       type: array
       items: File
     outputBinding:
       glob: '*.java'
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - tar
   - xf
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: untar


