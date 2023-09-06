:orphan:

HTSeq-Count
========================

``HTSeqCount`` · *1 contributor · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.htseq.htseqcount.versions import HTSeqCount_1_99_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "htseqcount_step",
           HTSeqCount_1_99_2(
               bams=None,
               gff_file=None,
           )
       )
       wf.output("out", source=htseqcount_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for HTSeqCount:

.. code-block:: bash

   # user inputs
   janis inputs HTSeqCount > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bams:
       - bams_0.bam
       - bams_1.bam
       gff_file: gff_file




5. Run HTSeqCount with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       HTSeqCount

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          HTSeqCount

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``HTSeqCount``
:URL: `https://htseq.readthedocs.io/en/release_0.11.1/count.html#count <https://htseq.readthedocs.io/en/release_0.11.1/count.html#count>`_
:Versions: 1.99.2
:Container: quay.io/biocontainers/htseq:1.99.2--py39haf81c86_0
:Authors: Jiaan Yu
:Citations: G Putri, S Anders, PT Pyl, JE Pimanda, F ZaniniAnalysing high-throughput sequencing data with HTSeq 2.0arXiv:2112.00939 (2021)
:DOI: https://htseq.readthedocs.io/en/release_0.11.1/overview.html
:Created: 2022-01-17
:Updated: 2022-01-17


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
======  ======  ===============


Additional configuration (inputs)
---------------------------------

========================  ==================  ===========================  ==========  ==================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                      type                prefix                         position  documentation
========================  ==================  ===========================  ==========  ==================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
bams                      Array<BAM>                                                3
gff_file                  File                                                      4
outputFilename            Optional<Filename>  >                                     5
format                    Optional<String>    --format=                             1  Format of the input data. Possible values are sam (for text SAM files) and bam (for binary BAM files). Default is sam.
order                     Optional<String>    --order=                              1  For paired-end data, the alignment have to be sorted either by read name or by alignment position. If your data is not sorted, use the samtools sort function of samtools to sort it. Use this option, with name or pos for <order> to indicate how the input data has been sorted. The default is name.If name is indicated, htseq-count expects all the alignments for the reads of a given read pair to appear in adjacent records in the input data. For pos, this is not expected; rather, read alignments whose mate alignment have not yet been seen are kept in a buffer in memory until the mate is found. While, strictly speaking, the latter will also work with unsorted data, sorting ensures that most alignment mates appear close to each other in the data and hence the buffer is much less likely to overflow.
max_reads_in_buffer       Optional<Integer>   --max-reads-in-buffer=                1  When <alignment_file> is paired end sorted by position, allow only so many reads to stay in memory until the mates are found (raising this number will use more memory). Has no effect for single end or paired end sorted by name. (default: 30000000)
stranded                  Optional<String>    --stranded=                           1  whether the data is from a strand-specific assay (default: yes)For stranded=no, a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature. For stranded=yes and single-end reads, the read has to be mapped to the same strand as the feature. For paired-end reads, the first read has to be on the same strand and the second read on the opposite strand. For stranded=reverse, these rules are reversed.
minaqual                  Optional<Integer>   -a                                    1  skip all reads with alignment quality lower than the given minimum value (default: 10 — Note: the default used to be 0 until version 0.5.4.)
type                      Optional<String>    --type=                               1  feature type (3rd column in GFF file) to be used, all features of other type are ignored (default, suitable for RNA-Seq analysis using an Ensembl GTF file: exon)
id                        Optional<String>    --idattr=                             1  GFF attribute to be used as feature ID. Several GFF lines with the same feature ID will be considered as parts of the same feature. The feature ID is used to identity the counts in the output table. The default, suitable for RNA-Seq analysis using an Ensembl GTF file, is gene_id.
additional_attr           Optional<String>    --additional-attr=                    1  Additional feature attributes, which will be printed as an additional column after the primary attribute column but before the counts column(s). The default is none, a suitable value to get gene names using an Ensembl GTF file is gene_name. To use more than one additional attribute, repeat the option in the command line more than once, with a single attribute each time, e.g. --additional-attr=gene_name --additional_attr=exon_number.
mode                      Optional<String>    --mode=                               1  Mode to handle reads overlapping more than one feature. Possible values for <mode> are union, intersection-strict and intersection-nonempty (default: union)
nonunique                 Optional<String>    --nonunique=                          1  Mode to handle reads that align to or are assigned to more than one feature in the overlap <mode> of choice (see -m option). <nonunique mode> are none and all (default: none)
secondary_alignments      Optional<String>    --secondary-alignments=               1  Mode to handle secondary alignments (SAM flag 0x100). <mode> can be score and ignore (default: score)
supplementary_alignments  Optional<String>    --supplementary-alignments=           1  Mode to handle supplementary/chimeric alignments (SAM flag 0x800). <mode> can be score and ignore (default: score)
samout                    Optional<String>    --samout=                             1  write out all SAM alignment records into an output SAM file called <samout>, annotating each line with its assignment to a feature or a special counter (as an optional field with tag ‘XF’)
========================  ==================  ===========================  ==========  ==================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task HTSeqCount {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       Array[File] bams
       File gff_file
       String? outputFilename
       String? format
       String? order
       Int? max_reads_in_buffer
       String? stranded
       Int? minaqual
       String? type
       String? id
       String? additional_attr
       String? mode
       String? nonunique
       String? secondary_alignments
       String? supplementary_alignments
       String? samout
     }

     command <<<
       set -e
       htseq-count \
         ~{if defined(format) then ("--format='" + format + "'") else ""} \
         ~{if defined(order) then ("--order='" + order + "'") else ""} \
         ~{if defined(max_reads_in_buffer) then ("--max-reads-in-buffer=" + max_reads_in_buffer) else ''} \
         ~{if defined(stranded) then ("--stranded='" + stranded + "'") else ""} \
         ~{if defined(minaqual) then ("-a " + minaqual) else ''} \
         ~{if defined(type) then ("--type='" + type + "'") else ""} \
         ~{if defined(id) then ("--idattr='" + id + "'") else ""} \
         ~{if defined(additional_attr) then ("--additional-attr='" + additional_attr + "'") else ""} \
         ~{if defined(mode) then ("--mode='" + mode + "'") else ""} \
         ~{if defined(nonunique) then ("--nonunique='" + nonunique + "'") else ""} \
         ~{if defined(secondary_alignments) then ("--secondary-alignments='" + secondary_alignments + "'") else ""} \
         ~{if defined(supplementary_alignments) then ("--supplementary-alignments='" + supplementary_alignments + "'") else ""} \
         ~{if defined(samout) then ("--samout='" + samout + "'") else ""} \
         ~{if length(bams) > 0 then "'" + sep("' '", bams) + "'" else ""} \
         '~{gff_file}' \
         > '~{select_first([outputFilename, "generated.htseq-count.txt"])}'
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "quay.io/biocontainers/htseq:1.99.2--py39haf81c86_0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }

     output {
       File out = select_first([outputFilename, "generated.htseq-count.txt"])
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: HTSeq-Count

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/htseq:1.99.2--py39haf81c86_0

   inputs:
   - id: bams
     label: bams
     type:
       type: array
       items: File
     inputBinding:
       position: 3
   - id: gff_file
     label: gff_file
     type: File
     inputBinding:
       position: 4
   - id: outputFilename
     label: outputFilename
     doc: ''
     type:
     - string
     - 'null'
     default: generated.htseq-count.txt
     inputBinding:
       prefix: '>'
       position: 5
   - id: format
     label: format
     doc: |-
       Format of the input data. Possible values are sam (for text SAM files) and bam (for binary BAM files). Default is sam.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --format=
       position: 1
       separate: false
   - id: order
     label: order
     doc: |-
       For paired-end data, the alignment have to be sorted either by read name or by alignment position. If your data is not sorted, use the samtools sort function of samtools to sort it. Use this option, with name or pos for <order> to indicate how the input data has been sorted. The default is name.If name is indicated, htseq-count expects all the alignments for the reads of a given read pair to appear in adjacent records in the input data. For pos, this is not expected; rather, read alignments whose mate alignment have not yet been seen are kept in a buffer in memory until the mate is found. While, strictly speaking, the latter will also work with unsorted data, sorting ensures that most alignment mates appear close to each other in the data and hence the buffer is much less likely to overflow.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --order=
       position: 1
       separate: false
   - id: max_reads_in_buffer
     label: max_reads_in_buffer
     doc: |-
       When <alignment_file> is paired end sorted by position, allow only so many reads to stay in memory until the mates are found (raising this number will use more memory). Has no effect for single end or paired end sorted by name. (default: 30000000)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-reads-in-buffer=
       position: 1
       separate: false
   - id: stranded
     label: stranded
     doc: |-
       whether the data is from a strand-specific assay (default: yes)For stranded=no, a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature. For stranded=yes and single-end reads, the read has to be mapped to the same strand as the feature. For paired-end reads, the first read has to be on the same strand and the second read on the opposite strand. For stranded=reverse, these rules are reversed.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --stranded=
       position: 1
       separate: false
   - id: minaqual
     label: minaqual
     doc: |-
       skip all reads with alignment quality lower than the given minimum value (default: 10 — Note: the default used to be 0 until version 0.5.4.)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -a
       position: 1
   - id: type
     label: type
     doc: |-
       feature type (3rd column in GFF file) to be used, all features of other type are ignored (default, suitable for RNA-Seq analysis using an Ensembl GTF file: exon)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --type=
       position: 1
       separate: false
   - id: id
     label: id
     doc: |-
       GFF attribute to be used as feature ID. Several GFF lines with the same feature ID will be considered as parts of the same feature. The feature ID is used to identity the counts in the output table. The default, suitable for RNA-Seq analysis using an Ensembl GTF file, is gene_id.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --idattr=
       position: 1
       separate: false
   - id: additional_attr
     label: additional_attr
     doc: |-
       Additional feature attributes, which will be printed as an additional column after the primary attribute column but before the counts column(s). The default is none, a suitable value to get gene names using an Ensembl GTF file is gene_name. To use more than one additional attribute, repeat the option in the command line more than once, with a single attribute each time, e.g. --additional-attr=gene_name --additional_attr=exon_number.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --additional-attr=
       position: 1
       separate: false
   - id: mode
     label: mode
     doc: |-
       Mode to handle reads overlapping more than one feature. Possible values for <mode> are union, intersection-strict and intersection-nonempty (default: union)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --mode=
       position: 1
       separate: false
   - id: nonunique
     label: nonunique
     doc: |-
       Mode to handle reads that align to or are assigned to more than one feature in the overlap <mode> of choice (see -m option). <nonunique mode> are none and all (default: none)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --nonunique=
       position: 1
       separate: false
   - id: secondary_alignments
     label: secondary_alignments
     doc: |-
       Mode to handle secondary alignments (SAM flag 0x100). <mode> can be score and ignore (default: score)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --secondary-alignments=
       position: 1
       separate: false
   - id: supplementary_alignments
     label: supplementary_alignments
     doc: |-
       Mode to handle supplementary/chimeric alignments (SAM flag 0x800). <mode> can be score and ignore (default: score)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --supplementary-alignments=
       position: 1
       separate: false
   - id: samout
     label: samout
     doc: |-
       write out all SAM alignment records into an output SAM file called <samout>, annotating each line with its assignment to a feature or a special counter (as an optional field with tag ‘XF’)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --samout=
       position: 1
       separate: false

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.htseq-count.txt
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - htseq-count
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: HTSeqCount


