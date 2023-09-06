:orphan:

BWA-Index
====================

``bwaIndex`` · *1 contributor · 1 version*

bwa - Burrows-Wheeler Alignment Tool
Index database sequences in the FASTA format.

Warning: `-a bwtsw' does not work for short genomes, while `-a is' and
         `-a div' do not work not for long genomes.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.bwa.index.versions import BwaIndex_0_7_15

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bwaindex_step",
           BwaIndex_0_7_15(
               reference=None,
           )
       )
       wf.output("out", source=bwaindex_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for bwaIndex:

.. code-block:: bash

   # user inputs
   janis inputs bwaIndex > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run bwaIndex with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bwaIndex

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          bwaIndex

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``bwaIndex``
:URL: `http://bio-bwa.sourceforge.net/bwa.shtml#3 <http://bio-bwa.sourceforge.net/bwa.shtml#3>`_
:Versions: v0.7.15
:Container: biocontainers/bwa:v0.7.15_cv3
:Authors: Michael Franklin
:Citations: The BWA-MEM algorithm has not been published yet.
:Created: 2020-02-14
:Updated: 2020-02-14


Outputs
-----------

======  ========  ===============
name    type      documentation
======  ========  ===============
out     FastaBwa
======  ========  ===============


Additional configuration (inputs)
---------------------------------

=========  =================  ========  ==========  =======================================================================
name       type               prefix      position  documentation
=========  =================  ========  ==========  =======================================================================
reference  Fasta                                 1
blockSize  Optional<Integer>  -b                    block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]
algorithm  Optional<String>   -a                    BWT construction algorithm: bwtsw, is or rb2 [auto]
                                                        - is	IS linear-time algorithm for constructing suffix array. It requires 5.37N memory where N is the size of the database. IS is moderately fast, but does not work with database larger than 2GB. IS is the default algorithm due to its simplicity. The current codes for IS algorithm are reimplemented by Yuta Mori.
                                                        - bwtsw	Algorithm implemented in BWT-SW. This method works with the whole human genome.
=========  =================  ========  ==========  =======================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task bwaIndex {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       File reference
       Int? blockSize
       String? algorithm
     }

     command <<<
       set -e
       cp -f '~{reference}' '.'
       bwa index \
         ~{if defined(blockSize) then ("-b " + blockSize) else ''} \
         ~{if defined(algorithm) then ("-a '" + algorithm + "'") else ""} \
         '~{basename(reference)}'
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "biocontainers/bwa:v0.7.15_cv3"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }

     output {
       File out = basename(reference)
       File out_amb = basename(reference) + ".amb"
       File out_ann = basename(reference) + ".ann"
       File out_bwt = basename(reference) + ".bwt"
       File out_pac = basename(reference) + ".pac"
       File out_sa = basename(reference) + ".sa"
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: BWA-Index

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: InitialWorkDirRequirement
     listing:
     - entry: $(inputs.reference)
   - class: DockerRequirement
     dockerPull: biocontainers/bwa:v0.7.15_cv3

   inputs:
   - id: reference
     label: reference
     type: File
     inputBinding:
       position: 1
   - id: blockSize
     label: blockSize
     doc: block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -b
   - id: algorithm
     label: algorithm
     doc: |
       BWT construction algorithm: bwtsw, is or rb2 [auto]
           - is	IS linear-time algorithm for constructing suffix array. It requires 5.37N memory where N is the size of the database. IS is moderately fast, but does not work with database larger than 2GB. IS is the default algorithm due to its simplicity. The current codes for IS algorithm are reimplemented by Yuta Mori.
           - bwtsw	Algorithm implemented in BWT-SW. This method works with the whole human genome.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -a

   outputs:
   - id: out
     label: out
     type: File
     secondaryFiles:
     - pattern: .amb
     - pattern: .ann
     - pattern: .bwt
     - pattern: .pac
     - pattern: .sa
     outputBinding:
       glob: $(inputs.reference.basename)
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - bwa
   - index
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: bwaIndex


