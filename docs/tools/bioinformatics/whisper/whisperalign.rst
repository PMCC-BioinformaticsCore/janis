:orphan:

Whisper-Align
============================

``whisperAlign`` · *1 contributor · 1 version*

Builds a whisper index


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.whisper.align.versions import WhisperAlign_2_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "whisperalign_step",
           WhisperAlign_2_0(
               index=None,
               fastq=None,
           )
       )
       wf.output("out", source=whisperalign_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for whisperAlign:

.. code-block:: bash

   # user inputs
   janis inputs whisperAlign > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fastq:
       - fastq_0.fastq.gz
       - fastq_1.fastq.gz
       index: index




5. Run whisperAlign with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       whisperAlign





Information
------------

:ID: ``whisperAlign``
:URL: `https://github.com/refresh-bio/Whisper <https://github.com/refresh-bio/Whisper>`_
:Versions: v2.0
:Container: drtomc/whisper:2.0
:Authors: Thomas Conway
:Citations: Deorowicz, S., Gudyś, A. (2019) Whisper 2: indel-sensitive short read mapping, biorXiv
:DOI: https://doi.org/10.1101/2019.12.18.881292
:Created: 2020-06-16
:Updated: 2020-06-16


Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     stdout<BAM>
======  ===========  ===============


Additional configuration (inputs)
---------------------------------

======  ===========  ========  ==========  ===========================
name    type         prefix      position  documentation
======  ===========  ========  ==========  ===========================
index   WhisperIdx                      2  base name for whisper index
fastq   FastqGzPair                     3  Paired end fastq reads
======  ===========  ========  ==========  ===========================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task whisperAlign {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File index
       File index_whisper_idxlut_long_dir
       File index_whisper_idxlut_long_rc
       File index_whisper_idxlut_short_dir
       File index_whisper_idxlut_short_rc
       File index_whisper_idxref_seq_desc
       File index_whisper_idxref_seq_dir_pck
       File index_whisper_idxref_seq_rc_pck
       File index_whisper_idxsa_dir
       File index_whisper_idxsa_rc
       Array[File] fastq
     }
     command <<<
       set -e
       whisper -stdout -t 4 -store-BAM \
         '~{index}' \
         ~{if length(fastq) > 0 then "'" + sep("' '", fastq) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "drtomc/whisper:2.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
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
   label: Whisper-Align
   doc: Builds a whisper index

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: drtomc/whisper:2.0

   inputs:
   - id: index
     label: index
     doc: base name for whisper index
     type: File
     secondaryFiles:
     - pattern: .whisper_idx.lut_long_dir
     - pattern: .whisper_idx.lut_long_rc
     - pattern: .whisper_idx.lut_short_dir
     - pattern: .whisper_idx.lut_short_rc
     - pattern: .whisper_idx.ref_seq_desc
     - pattern: .whisper_idx.ref_seq_dir_pck
     - pattern: .whisper_idx.ref_seq_rc_pck
     - pattern: .whisper_idx.sa_dir
     - pattern: .whisper_idx.sa_rc
     inputBinding:
       position: 2
   - id: fastq
     label: fastq
     doc: Paired end fastq reads
     type:
       type: array
       items: File
     inputBinding:
       position: 3

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - whisper
   - -stdout
   - -t
   - '4'
   - -store-BAM
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: whisperAlign


