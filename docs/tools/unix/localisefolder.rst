:orphan:

LocaliseFolder
==============

``LocaliseFolder`` · *0 contributors · 1 version*

localise folder


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.localisefolder import LocaliseFolder

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "localisefolder_step",
           LocaliseFolder(
               dir=None,
           )
       )
       wf.output("out", source=localisefolder_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for LocaliseFolder:

.. code-block:: bash

   # user inputs
   janis inputs LocaliseFolder > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       dir: null




5. Run LocaliseFolder with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       LocaliseFolder

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          LocaliseFolder

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``LocaliseFolder``
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

======  =========  ========  ==========  ===============
name    type       prefix      position  documentation
======  =========  ========  ==========  ===============
dir     Directory                     3
======  =========  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task LocaliseFolder {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       Directory dir
     }

     command <<<
       set -e
        \
         cp \
         -r \
         '~{dir}' \
         .
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
       Directory out = dir
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: LocaliseFolder

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

   inputs:
   - id: dir
     label: dir
     type: Directory
     inputBinding:
       position: 3

   outputs:
   - id: out
     label: out
     type: Directory
     outputBinding:
       glob: $(inputs.dir)
       loadContents: false
   stdout: _stdout
   stderr: _stderr
   arguments:
   - position: 1
     valueFrom: cp
     shellQuote: false
   - position: 2
     valueFrom: -r
     shellQuote: false
   - position: 4
     valueFrom: .
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: LocaliseFolder


