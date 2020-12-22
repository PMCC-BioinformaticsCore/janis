:orphan:

SamTools: Index
===============================

``SamToolsIndex`` · *1 contributor · 2 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.samtools.index.versions import SamToolsIndex_1_9

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "samtoolsindex_step",
           SamToolsIndex_1_9(
               bam=None,
           )
       )
       wf.output("out", source=samtoolsindex_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SamToolsIndex:

.. code-block:: bash

   # user inputs
   janis inputs SamToolsIndex > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: null




5. Run SamToolsIndex with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SamToolsIndex





Information
------------

:ID: ``SamToolsIndex``
:URL: `http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS <http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS>`_
:Versions: 1.9.0, 1.7.0
:Container: quay.io/biocontainers/samtools:1.9--h8571acd_11
:Authors: Michael Franklin
:Citations: None
:Created: 2019-12-17
:Updated: 2019-12-17


Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedBam
======  ==========  ===============


Additional configuration (inputs)
---------------------------------

=======  =====================  ========  ==========  ===============
name     type                   prefix      position  documentation
=======  =====================  ========  ==========  ===============
bam      Union<BAM, SAM, CRAM>                    10
threads  Optional<Integer>      -@                10
=======  =====================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task SamToolsIndex {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File bam
       Int? threads
     }
     command <<<
       set -e
       samtools index \
         '-b' \
         ~{basename(bam)} \
         ~{if defined(select_first([threads, select_first([runtime_cpu, 1])])) then ("-@ " + select_first([threads, select_first([runtime_cpu, 1])])) else ''}
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
       File out = basename(bam)
       File out_bai = basename(bam) + ".bai"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'SamTools: Index'
   doc: ''

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/samtools:1.9--h8571acd_11

   inputs:
   - id: bam
     label: bam
     type: File
     inputBinding:
       position: 10
   - id: threads
     label: threads
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -@
       position: 10
       valueFrom: $([inputs.runtime_cpu, 1].filter(function (inner) { return inner !=
         null })[0])

   outputs:
   - id: out
     label: out
     type: File
     secondaryFiles:
     - pattern: .bai
     outputBinding:
       glob: $(inputs.bam)
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - samtools
   - index
   arguments:
   - position: 4
     valueFrom: -b

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: SamToolsIndex


