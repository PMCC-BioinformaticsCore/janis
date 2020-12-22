:orphan:

Sequenza: bam2seqz
=================================

``SeqzBam2Seqz`` · *2 contributors · 2 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.sequenza.bam2seqz.versions import SequenzaBam2Seqz_2_2_0_9000

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "seqzbam2seqz_step",
           SequenzaBam2Seqz_2_2_0_9000(
               normal=None,
               tumour=None,
               wiggle_file=None,
               fasta_reference=None,
           )
       )
       wf.output("out", source=seqzbam2seqz_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SeqzBam2Seqz:

.. code-block:: bash

   # user inputs
   janis inputs SeqzBam2Seqz > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fasta_reference: fasta_reference.fasta
       normal: normal.bam
       tumour: tumour.bam
       wiggle_file: wiggle_file




5. Run SeqzBam2Seqz with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SeqzBam2Seqz





Information
------------

:ID: ``SeqzBam2Seqz``
:URL: `http://www.cbs.dtu.dk/biotools/sequenza/ <http://www.cbs.dtu.dk/biotools/sequenza/>`_
:Versions: 3.0.0, 2.2.0.9000
:Container: sequenza/sequenza:2.2.0.9000
:Authors: mumbler, evanwehi
:Citations: None
:Created: 2019-12-10
:Updated: 2019-12-10


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
======  ======  ===============


Additional configuration (inputs)
---------------------------------

===============  ==================  ========  ==========  ==============================================================================================
name             type                prefix      position  documentation
===============  ==================  ========  ==========  ==============================================================================================
normal           IndexedBam          --normal           2  Name of the BAM/pileup file from the reference/normal sample
tumour           IndexedBam          --tumor            4  Name of the BAM/pileup file from the reference/normal sample
wiggle_file      File                -gc                6  The GC-content wiggle file
fasta_reference  FastaFai            --fasta            8  The reference FASTA file used to generate the intermediate pileup. Required when input are BAM
output_filename  Optional<Filename>  --output          10  Name of the output file. To use gzip compression name the file ending in .gz. Default STDOUT.
===============  ==================  ========  ==========  ==============================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task SeqzBam2Seqz {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File normal
       File normal_bai
       File tumour
       File tumour_bai
       File wiggle_file
       File fasta_reference
       File fasta_reference_fai
       String? output_filename
     }
     command <<<
       set -e
       sequenza-utils bam2seqz \
         --normal '~{normal}' \
         --tumor '~{tumour}' \
         -gc '~{wiggle_file}' \
         --fasta '~{fasta_reference}' \
         --output '~{select_first([output_filename, "generated.gz"])}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "sequenza/sequenza:2.2.0.9000"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([output_filename, "generated.gz"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'Sequenza: bam2seqz'
   doc: ''

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: sequenza/sequenza:2.2.0.9000

   inputs:
   - id: normal
     label: normal
     doc: Name of the BAM/pileup file from the reference/normal sample
     type: File
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       prefix: --normal
       position: 2
   - id: tumour
     label: tumour
     doc: Name of the BAM/pileup file from the reference/normal sample
     type: File
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       prefix: --tumor
       position: 4
   - id: wiggle_file
     label: wiggle_file
     doc: The GC-content wiggle file
     type: File
     inputBinding:
       prefix: -gc
       position: 6
   - id: fasta_reference
     label: fasta_reference
     doc: |-
       The reference FASTA file used to generate the intermediate pileup. Required when input are BAM
     type: File
     secondaryFiles:
     - pattern: .fai
     inputBinding:
       prefix: --fasta
       position: 8
   - id: output_filename
     label: output_filename
     doc: |-
       Name of the output file. To use gzip compression name the file ending in .gz. Default STDOUT.
     type:
     - string
     - 'null'
     default: generated.gz
     inputBinding:
       prefix: --output
       position: 10

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.gz
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - sequenza-utils
   - bam2seqz
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: SeqzBam2Seqz


