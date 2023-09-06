:orphan:

LocaliseFastaWithDict
=====================

``LocaliseFastaWithDict`` · *1 contributor · 1 version*

Localise FastaWithDict
        


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.localisefastawithdict import LocaliseFastaWithDict

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "localisefastawithdict_step",
           LocaliseFastaWithDict(
               reference=None,
           )
       )
       wf.output("out", source=localisefastawithdict_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for LocaliseFastaWithDict:

.. code-block:: bash

   # user inputs
   janis inputs LocaliseFastaWithDict > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run LocaliseFastaWithDict with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       LocaliseFastaWithDict

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          LocaliseFastaWithDict

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``LocaliseFastaWithDict``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Container: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956
:Authors: Jiaan Yu
:Citations: None
:Created: 2022-01-07 00:00:00
:Updated: 2022-01-07 00:00:00


Outputs
-----------

======  ================  ===============
name    type              documentation
======  ================  ===============
out     FastaWithIndexes
======  ================  ===============


Additional configuration (inputs)
---------------------------------

=========  ================  ========  ==========  ===============
name       type              prefix    position    documentation
=========  ================  ========  ==========  ===============
reference  FastaWithIndexes
=========  ================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task LocaliseFastaWithDict {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
     }

     command <<<
       set -e
       cp -f '~{reference}' '.'
       cp -f '~{reference_fai}' .
       cp -f '~{reference_amb}' .
       cp -f '~{reference_ann}' .
       cp -f '~{reference_bwt}' .
       cp -f '~{reference_pac}' .
       cp -f '~{reference_sa}' .
       cp -f '~{reference_dict}' .
    
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
       File out = basename(basename(reference))
       File out_fai = basename(basename(reference)) + ".fai"
       File out_amb = basename(basename(reference)) + ".amb"
       File out_ann = basename(basename(reference)) + ".ann"
       File out_bwt = basename(basename(reference)) + ".bwt"
       File out_pac = basename(basename(reference)) + ".pac"
       File out_sa = basename(basename(reference)) + ".sa"
       File out_dict = sub(sub(sub(basename(basename(reference)), "\\.fasta$", ".dict"), "\\.fna$", ".dict"), "\\.fa$", ".dict")
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: LocaliseFastaWithDict

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: InitialWorkDirRequirement
     listing:
     - entry: $(inputs.reference)
   - class: DockerRequirement
     dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

   inputs:
   - id: reference
     label: reference
     type: File
     secondaryFiles:
     - pattern: .fai
     - pattern: .amb
     - pattern: .ann
     - pattern: .bwt
     - pattern: .pac
     - pattern: .sa
     - pattern: ^.dict

   outputs:
   - id: out
     label: out
     type: File
     secondaryFiles:
     - pattern: .fai
     - pattern: .amb
     - pattern: .ann
     - pattern: .bwt
     - pattern: .pac
     - pattern: .sa
     - pattern: ^.dict
     outputBinding:
       glob: $(inputs.reference.basename.basename)
       loadContents: false
   stdout: _stdout
   stderr: _stderr
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: LocaliseFastaWithDict


