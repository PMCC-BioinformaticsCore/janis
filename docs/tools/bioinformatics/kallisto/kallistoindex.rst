:orphan:

Kallisto-Index
==============================

``kallistoIndex`` · *1 contributor · 1 version*

Builds a kallisto index


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.kallisto.index.versions import KallistoIndex_0_46_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "kallistoindex_step",
           KallistoIndex_0_46_2(
               reference=None,
           )
       )
       wf.output("out", source=kallistoindex_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for kallistoIndex:

.. code-block:: bash

   # user inputs
   janis inputs kallistoIndex > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run kallistoIndex with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       kallistoIndex





Information
------------

:ID: ``kallistoIndex``
:URL: `https://pachterlab.github.io/kallisto/manual.html <https://pachterlab.github.io/kallisto/manual.html>`_
:Versions: v0.46.2
:Container: quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1
:Authors: Thomas Conway
:Citations: NL Bray, H Pimentel, P Melsted and L Pachter, Near optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, p 525--527 (2016).
:DOI: https://doi.org/10.1038/nbt.3519
:Created: 2020-05-25
:Updated: 2020-05-25


Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     KallistoIdx
======  ===========  ===============


Additional configuration (inputs)
---------------------------------

=========  ==================  ========  ==========  =================================================
name       type                prefix      position  documentation
=========  ==================  ========  ==========  =================================================
reference  Fasta                                  3  Filename for a reference transcriptome
kmer_size  Optional<Integer>   -k                 1  k-mer (odd) length (default: 31, max value: 31)
index      Optional<Filename>  -i                 2  Filename for the kallisto index to be constructed
=========  ==================  ========  ==========  =================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task kallistoIndex {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Int? kmer_size
       String? index
       File reference
     }
     command <<<
       set -e
       cp -f '~{reference}' '.'
       kallisto index \
         ~{if defined(kmer_size) then ("-k " + kmer_size) else ''} \
         -i '~{select_first([index, "generated.kidx"])}' \
         '~{basename(reference)}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 2, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([index, "generated.kidx"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Kallisto-Index
   doc: Builds a kallisto index

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: InitialWorkDirRequirement
     listing:
     - entry: $(inputs.reference)
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1

   inputs:
   - id: kmer_size
     label: kmer_size
     doc: 'k-mer (odd) length (default: 31, max value: 31)'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -k
       position: 1
   - id: index
     label: index
     doc: Filename for the kallisto index to be constructed
     type:
     - string
     - 'null'
     default: generated.kidx
     inputBinding:
       prefix: -i
       position: 2
   - id: reference
     label: reference
     doc: Filename for a reference transcriptome
     type: File
     inputBinding:
       position: 3

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.kidx
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - kallisto
   - index
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: kallistoIndex


