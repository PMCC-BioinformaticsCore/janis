:orphan:

LocaliseFastqGzPair
===================

``LocaliseFastqGzPair`` · *1 contributor · 1 version*

Localise Arrary of FastqGZ pairs
        


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.localisefastqgzpair import LocaliseFastqGzPair

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "localisefastqgzpair_step",
           LocaliseFastqGzPair(
               fastqs=None,
           )
       )
       wf.output("out", source=localisefastqgzpair_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for LocaliseFastqGzPair:

.. code-block:: bash

   # user inputs
   janis inputs LocaliseFastqGzPair > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fastqs:
       - fastqs_0.fastq.gz
       - fastqs_1.fastq.gz




5. Run LocaliseFastqGzPair with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       LocaliseFastqGzPair

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          LocaliseFastqGzPair

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``LocaliseFastqGzPair``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Container: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956
:Authors: Jiaan Yu
:Citations: None
:Created: 2022-01-07 00:00:00
:Updated: 2022-01-07 00:00:00


Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     FastqGzPair
======  ===========  ===============


Additional configuration (inputs)
---------------------------------

===========  ==================  ========  ==========  ===============
name         type                prefix      position  documentation
===========  ==================  ========  ==========  ===============
fastqs       FastqGzPair                            2
outputRead1  Optional<Filename>
outputRead2  Optional<Filename>
===========  ==================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task LocaliseFastqGzPair {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       Array[File] fastqs
       String? outputRead1
       String? outputRead2
     }

     command <<<
       set -e
        \
         cp \
         ~{if length(fastqs) > 0 then "'" + sep("' '", fastqs) + "'" else ""} \
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
       Array[File] out = [basename(select_first([outputRead1, "~{fastqs[0]}"])), basename(select_first([outputRead2, "~{fastqs[1]}"]))]
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: LocaliseFastqGzPair

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

   inputs:
   - id: fastqs
     label: fastqs
     type:
       type: array
       items: File
     inputBinding:
       position: 2
   - id: outputRead1
     label: outputRead1
     type:
     - string
     - 'null'
     default: generated
   - id: outputRead2
     label: outputRead2
     type:
     - string
     - 'null'
     default: generated

   outputs:
   - id: out
     label: out
     type:
       type: array
       items: File
     outputBinding:
       glob:
       - $(inputs.outputRead1.basename)
       - $(inputs.outputRead2.basename)
       loadContents: false
   stdout: _stdout
   stderr: _stderr
   arguments:
   - position: 1
     valueFrom: cp
     shellQuote: false
   - position: 3
     valueFrom: .
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: LocaliseFastqGzPair


