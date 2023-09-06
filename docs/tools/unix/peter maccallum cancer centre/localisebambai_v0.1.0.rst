:orphan:

LocaliseBamBai
==============

``LocaliseBamBai`` · *1 contributor · 1 version*

Localise BamBai
        


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.localisebambai import LocaliseBamBai

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "localisebambai_step",
           LocaliseBamBai(
               bam=None,
           )
       )
       wf.output("out", source=localisebambai_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for LocaliseBamBai:

.. code-block:: bash

   # user inputs
   janis inputs LocaliseBamBai > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam




5. Run LocaliseBamBai with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       LocaliseBamBai

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          LocaliseBamBai

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``LocaliseBamBai``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Container: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956
:Authors: Jiaan Yu
:Citations: None
:Created: 2022-01-07 00:00:00
:Updated: 2022-01-07 00:00:00


Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedBam
======  ==========  ===============


Additional configuration (inputs)
---------------------------------

======  ==========  ========  ==========  ===============
name    type        prefix    position    documentation
======  ==========  ========  ==========  ===============
bam     IndexedBam
======  ==========  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task LocaliseBamBai {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       File bam
       File bam_bai
     }

     command <<<
       set -e
       cp -f '~{bam}' '.'
       cp -f '~{bam_bai}' .
    
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
       File out = basename(basename(bam))
       File out_bai = basename(basename(bam)) + ".bai"
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: LocaliseBamBai

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: InitialWorkDirRequirement
     listing:
     - entry: $(inputs.bam)
   - class: DockerRequirement
     dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

   inputs:
   - id: bam
     label: bam
     type: File
     secondaryFiles:
     - pattern: .bai

   outputs:
   - id: out
     label: out
     type: File
     secondaryFiles:
     - pattern: .bai
     outputBinding:
       glob: $(inputs.bam.basename.basename)
       loadContents: false
   stdout: _stdout
   stderr: _stderr
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: LocaliseBamBai


