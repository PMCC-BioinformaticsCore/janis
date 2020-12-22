:orphan:

Whisper-Index
============================

``WhisperIndex`` · *1 contributor · 1 version*

Builds a whisper index


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.whisper.index.versions import WhisperIndex_2_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "whisperindex_step",
           WhisperIndex_2_0(
               index_name=None,
               fasta=None,
           )
       )
       wf.output("out", source=whisperindex_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for WhisperIndex:

.. code-block:: bash

   # user inputs
   janis inputs WhisperIndex > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fasta:
       - fasta_0.fasta
       - fasta_1.fasta
       index_name: <value>




5. Run WhisperIndex with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       WhisperIndex





Information
------------

:ID: ``WhisperIndex``
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

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     WhisperIdx
======  ==========  ===============


Additional configuration (inputs)
---------------------------------

==========  ============  ========  ==========  ====================
name        type          prefix      position  documentation
==========  ============  ========  ==========  ====================
index_name  String                           2  name of the index
fasta       Array<Fasta>                     3  FASTA files to index
==========  ============  ========  ==========  ====================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task WhisperIndex {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       String index_name
       Array[File] fasta
     }
     command <<<
       set -e
       whindex whisper-index \
         '~{index_name}' \
         ~{if length(fasta) > 0 then "'" + sep("' '", fasta) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "drtomc/whisper:2.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = ("whisper-index/" + index_name)
       File out_whisper_idxlut_long_dir = ("whisper-index/" + index_name) + ".whisper_idx.lut_long_dir"
       File out_whisper_idxlut_long_rc = ("whisper-index/" + index_name) + ".whisper_idx.lut_long_rc"
       File out_whisper_idxlut_short_dir = ("whisper-index/" + index_name) + ".whisper_idx.lut_short_dir"
       File out_whisper_idxlut_short_rc = ("whisper-index/" + index_name) + ".whisper_idx.lut_short_rc"
       File out_whisper_idxref_seq_desc = ("whisper-index/" + index_name) + ".whisper_idx.ref_seq_desc"
       File out_whisper_idxref_seq_dir_pck = ("whisper-index/" + index_name) + ".whisper_idx.ref_seq_dir_pck"
       File out_whisper_idxref_seq_rc_pck = ("whisper-index/" + index_name) + ".whisper_idx.ref_seq_rc_pck"
       File out_whisper_idxsa_dir = ("whisper-index/" + index_name) + ".whisper_idx.sa_dir"
       File out_whisper_idxsa_rc = ("whisper-index/" + index_name) + ".whisper_idx.sa_rc"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Whisper-Index
   doc: Builds a whisper index

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: drtomc/whisper:2.0

   inputs:
   - id: index_name
     label: index_name
     doc: name of the index
     type: string
     inputBinding:
       position: 2
   - id: fasta
     label: fasta
     doc: FASTA files to index
     type:
       type: array
       items: File
     inputBinding:
       position: 3

   outputs:
   - id: out
     label: out
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
     outputBinding:
       glob: $(("whisper-index/" + inputs.index_name))
       outputEval: $(("whisper-index/" + inputs.index_name.basename))
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - whindex
   - whisper-index
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: WhisperIndex


