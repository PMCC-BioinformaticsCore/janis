:orphan:

SamTools: faidx
===============================

``SamToolsFaidx`` · *1 contributor · 2 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.samtools.faidx.versions import SamToolsFaidx_1_9

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "samtoolsfaidx_step",
           SamToolsFaidx_1_9(
               reference=None,
           )
       )
       wf.output("out", source=samtoolsfaidx_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SamToolsFaidx:

.. code-block:: bash

   # user inputs
   janis inputs SamToolsFaidx > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run SamToolsFaidx with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SamToolsFaidx





Information
------------

:ID: ``SamToolsFaidx``
:URL: `http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS <http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS>`_
:Versions: 1.9.0, 1.7.0
:Container: quay.io/biocontainers/samtools:1.9--h8571acd_11
:Authors: Michael Franklin
:Citations: None
:Created: 2020-02-14
:Updated: 2020-02-14


Outputs
-----------

======  ========  ===============
name    type      documentation
======  ========  ===============
out     FastaFai
======  ========  ===============


Additional configuration (inputs)
---------------------------------

=========  ======  ========  ==========  ===============
name       type    prefix      position  documentation
=========  ======  ========  ==========  ===============
reference  Fasta                      1
=========  ======  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task SamToolsFaidx {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File reference
     }
     command <<<
       set -e
       cp -f '~{reference}' '.'
       samtools faidx \
         '~{basename(reference)}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/samtools:1.9--h8571acd_11"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = basename(reference)
       File out_fai = basename(reference) + ".fai"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'SamTools: faidx'
   doc: ''

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: InitialWorkDirRequirement
     listing:
     - entry: $(inputs.reference)
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/samtools:1.9--h8571acd_11

   inputs:
   - id: reference
     label: reference
     type: File
     inputBinding:
       position: 1

   outputs:
   - id: out
     label: out
     type: File
     secondaryFiles:
     - pattern: .fai
     outputBinding:
       glob: $(inputs.reference.basename)
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - samtools
   - faidx
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: SamToolsFaidx


