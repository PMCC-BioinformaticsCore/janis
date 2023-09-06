:orphan:

RNASeqQC
========

``RNASeqQC`` · *1 contributor · 1 version*

Usage: rnaseqc [gtf] [bam] [output] \{OPTIONS\}



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.rnaseqqc.versions import RNASeqQC_2_3_5

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "rnaseqqc_step",
           RNASeqQC_2_3_5(
               gtf=None,
               bam=None,
           )
       )
       wf.output("out_gene_fragments", source=rnaseqqc_step.out_gene_fragments)
       wf.output("out_gene_reads", source=rnaseqqc_step.out_gene_reads)
       wf.output("out_gene_tpm", source=rnaseqqc_step.out_gene_tpm)
       wf.output("out_metrics_tsv", source=rnaseqqc_step.out_metrics_tsv)
       wf.output("out_coverage_tsv", source=rnaseqqc_step.out_coverage_tsv)
       wf.output("out_exon_reads", source=rnaseqqc_step.out_exon_reads)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for RNASeqQC:

.. code-block:: bash

   # user inputs
   janis inputs RNASeqQC > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       gtf: gtf




5. Run RNASeqQC with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       RNASeqQC

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          RNASeqQC

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``RNASeqQC``
:URL: `https://github.com/getzlab/rnaseqc <https://github.com/getzlab/rnaseqc>`_
:Versions: 2.3.5
:Container: quay.io/biocontainers/rna-seqc:2.3.5--he24ac62_2
:Authors: Jiaan Yu
:Citations: None
:Created: 2021-09-10
:Updated: 2021-10-19


Outputs
-----------

==================  =============  ===============
name                type           documentation
==================  =============  ===============
out_gene_fragments  File
out_gene_reads      File
out_gene_tpm        File
out_metrics_tsv     tsv
out_coverage_tsv    Optional<tsv>
out_exon_reads      File
==================  =============  ===============


Additional configuration (inputs)
---------------------------------

===================  =================  =====================  ==========  ========================================================================================================================================================================================================
name                 type               prefix                   position  documentation
===================  =================  =====================  ==========  ========================================================================================================================================================================================================
gtf                  File                                               1  The input GTF file containing features to check the bam against
bam                  IndexedBam                                         2  The input SAM/BAM file containing reads to process
output_dir           Optional<String>                                   3  Output directory
sample               Optional<String>   --sample                        4  The name of the current sample. Default: The bam's filename
bed                  Optional<bed>      --bed                           4  Optional input BED file containing non-overlapping exons used for fragment size calculations
fasta                Optional<Fasta>    --fasta                         4  Optional input FASTA/FASTQ file containing the reference sequence used for parsing CRAM files
chimeric_distance    Optional<Integer>  --chimeric-distance             4  Set the maximum accepted distance between read mates. Mates beyond this distance will be counted as chimeric pairs. Default: 2000000 [bp]
fragment_samples     Optional<Integer>  --fragment-samples              4  Set the number of samples to take when computing fragment sizes. Requires the --bed argument. Default: 1000000
mapping_quality      Optional<Integer>  --mapping-quality               4  Set the lower bound on read quality for exon coverage counting. Reads below this number are excluded from coverage metrics. Default: 255
base_mismatch        Optional<Integer>  --base-mismatch                 4  Set the maximum number of allowed mismatches between a read and the reference sequence. Reads with more than this number of mismatches are excluded from coverage metrics. Default: 6
offset               Optional<Integer>  --offset                        4  Set the offset into the gene for the 3' and 5' windows in bias calculation. A positive value shifts the 3' and 5' windows towards eachother, while a negative value shifts them apart. Default: 150 [bp]
window_size          Optional<Integer>  --window-size                   4  Set the size of the 3' and 5' windows in bias calculation. Default: 100 [bp]
gene_length          Optional<Integer>  --gene-length                   4  Set the minimum size of a gene for bias calculation. Genes below this size are ignored in the calculation. Default: 600 [bp]
legacy               Optional<Boolean>  --legacy                        4  Use legacy counting rules. Gene and exon counts match output of RNA-SeQC 1.1.9
stranded             Optional<String>   --stranded                      4  Use strand-specific metrics. Only features on the same strand of a read will be considered. Allowed values are 'RF', 'rf', 'FR', and 'fr'
verbose              Optional<Boolean>  --verbose                       4  Give some feedback about what's going on. Supply this argument twice for progress updates while parsing the bam
tag                  Optional<String>   --tag                           4  Filter out reads with the specified tag.
chimeric_tag         Optional<String>   --chimeric-tag                  4  Reads maked with the specified tag will be labeled as Chimeric. Defaults to 'mC' for STAR
exclude_chimeric     Optional<Boolean>  --exclude-chimeric              4  Exclude chimeric reads from the read counts
unpaired             Optional<Boolean>  --unpaired                      4  Allow unpaired reads to be quantified. Required for single-end libraries
rpkm                 Optional<Boolean>  --rpkm                          4  Output gene RPKM values instead of TPMs
coverage             Optional<Boolean>  --coverage                      4  If this flag is provided, coverage statistics for each transcript will be written to a table. Otherwise, only summary coverage statistics are generated and added to the metrics table
coverage_mask        Optional<Integer>  --coverage-mask                 4  Sets how many bases at both ends of a transcript are masked out when computing per-base exon coverage. Default: 500bp
detection_threshold  Optional<Integer>  --detection-threshold           4  Number of counts on a gene to consider the gene 'detected'. Additionally, genes below this limit are excluded from 3' bias computation. Default: 5 reads
===================  =================  =====================  ==========  ========================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task RNASeqQC {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       File gtf
       File bam
       File bam_bai
       String? output_dir
       String? sample
       File? bed
       File? fasta
       Int? chimeric_distance
       Int? fragment_samples
       Int? mapping_quality
       Int? base_mismatch
       Int? offset
       Int? window_size
       Int? gene_length
       Boolean? legacy
       String? stranded
       Boolean? verbose
       String? tag
       String? chimeric_tag
       Boolean? exclude_chimeric
       Boolean? unpaired
       Boolean? rpkm
       Boolean? coverage
       Int? coverage_mask
       Int? detection_threshold
     }

     command <<<
       set -e
       rnaseqc \
         '~{gtf}' \
         '~{bam}' \
         ~{if defined(select_first([output_dir, "."])) then ("'" + select_first([output_dir, "."]) + "'") else ""} \
         ~{if defined(sample) then ("--sample '" + sample + "'") else ""} \
         ~{if defined(bed) then ("--bed '" + bed + "'") else ""} \
         ~{if defined(fasta) then ("--fasta '" + fasta + "'") else ""} \
         ~{if defined(chimeric_distance) then ("--chimeric-distance " + chimeric_distance) else ''} \
         ~{if defined(fragment_samples) then ("--fragment-samples " + fragment_samples) else ''} \
         ~{if defined(mapping_quality) then ("--mapping-quality " + mapping_quality) else ''} \
         ~{if defined(base_mismatch) then ("--base-mismatch " + base_mismatch) else ''} \
         ~{if defined(offset) then ("--offset " + offset) else ''} \
         ~{if defined(window_size) then ("--window-size " + window_size) else ''} \
         ~{if defined(gene_length) then ("--gene-length " + gene_length) else ''} \
         ~{if (defined(legacy) && select_first([legacy])) then "--legacy" else ""} \
         ~{if defined(stranded) then ("--stranded '" + stranded + "'") else ""} \
         ~{if (defined(verbose) && select_first([verbose])) then "--verbose" else ""} \
         ~{if defined(tag) then ("--tag '" + tag + "'") else ""} \
         ~{if defined(chimeric_tag) then ("--chimeric-tag '" + chimeric_tag + "'") else ""} \
         ~{if (defined(exclude_chimeric) && select_first([exclude_chimeric])) then "--exclude-chimeric" else ""} \
         ~{if (defined(unpaired) && select_first([unpaired])) then "--unpaired" else ""} \
         ~{if (defined(rpkm) && select_first([rpkm])) then "--rpkm" else ""} \
         ~{if (defined(coverage) && select_first([coverage])) then "--coverage" else ""} \
         ~{if defined(coverage_mask) then ("--coverage-mask " + coverage_mask) else ''} \
         ~{if defined(detection_threshold) then ("--detection-threshold " + detection_threshold) else ''}
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "quay.io/biocontainers/rna-seqc:2.3.5--he24ac62_2"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4, 4])}G"
       preemptible: 2
     }

     output {
       File out_gene_fragments = "~{select_first([output_dir, "."])}/~{sample}.gene_fragments.gct"
       File out_gene_reads = "~{select_first([output_dir, "."])}/~{sample}.gene_reads.gct"
       File out_gene_tpm = "~{select_first([output_dir, "."])}/~{sample}.gene_tpm.gct"
       File out_metrics_tsv = "~{select_first([output_dir, "."])}/~{sample}.metrics.tsv"
       File? out_coverage_tsv = "~{select_first([output_dir, "."])}/~{sample}.coverage.tsv"
       File out_exon_reads = "~{select_first([output_dir, "."])}/~{sample}.exon_reads.gct"
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: RNASeqQC

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/rna-seqc:2.3.5--he24ac62_2

   inputs:
   - id: gtf
     label: gtf
     doc: The input GTF file containing features to check the bam against
     type: File
     inputBinding:
       position: 1
   - id: bam
     label: bam
     doc: The input SAM/BAM file containing reads to process
     type: File
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       position: 2
   - id: output_dir
     label: output_dir
     doc: Output directory
     type: string
     default: .
     inputBinding:
       position: 3
   - id: sample
     label: sample
     doc: "The name of the current sample. Default: The bam's filename"
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sample
       position: 4
   - id: bed
     label: bed
     doc: |-
       Optional input BED file containing non-overlapping exons used for fragment size calculations
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --bed
       position: 4
   - id: fasta
     label: fasta
     doc: |-
       Optional input FASTA/FASTQ file containing the reference sequence used for parsing CRAM files
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --fasta
       position: 4
   - id: chimeric_distance
     label: chimeric_distance
     doc: |-
       Set the maximum accepted distance between read mates. Mates beyond this distance will be counted as chimeric pairs. Default: 2000000 [bp]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --chimeric-distance
       position: 4
   - id: fragment_samples
     label: fragment_samples
     doc: |-
       Set the number of samples to take when computing fragment sizes. Requires the --bed argument. Default: 1000000
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --fragment-samples
       position: 4
   - id: mapping_quality
     label: mapping_quality
     doc: |-
       Set the lower bound on read quality for exon coverage counting. Reads below this number are excluded from coverage metrics. Default: 255
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --mapping-quality
       position: 4
   - id: base_mismatch
     label: base_mismatch
     doc: |-
       Set the maximum number of allowed mismatches between a read and the reference sequence. Reads with more than this number of mismatches are excluded from coverage metrics. Default: 6
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --base-mismatch
       position: 4
   - id: offset
     label: offset
     doc: |2-
        Set the offset into the gene for the 3' and 5' windows in bias calculation. A positive value shifts the 3' and 5' windows towards eachother, while a negative value shifts them apart. Default: 150 [bp]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --offset
       position: 4
   - id: window_size
     label: window_size
     doc: "Set the size of the 3' and 5' windows in bias calculation. Default: 100 [bp]"
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --window-size
       position: 4
   - id: gene_length
     label: gene_length
     doc: |-
       Set the minimum size of a gene for bias calculation. Genes below this size are ignored in the calculation. Default: 600 [bp]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --gene-length
       position: 4
   - id: legacy
     label: legacy
     doc: Use legacy counting rules. Gene and exon counts match output of RNA-SeQC 1.1.9
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --legacy
       position: 4
   - id: stranded
     label: stranded
     doc: |-
       Use strand-specific metrics. Only features on the same strand of a read will be considered. Allowed values are 'RF', 'rf', 'FR', and 'fr'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --stranded
       position: 4
   - id: verbose
     label: verbose
     doc: |-
       Give some feedback about what's going on. Supply this argument twice for progress updates while parsing the bam
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --verbose
       position: 4
   - id: tag
     label: tag
     doc: Filter out reads with the specified tag.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --tag
       position: 4
   - id: chimeric_tag
     label: chimeric_tag
     doc: |-
       Reads maked with the specified tag will be labeled as Chimeric. Defaults to 'mC' for STAR
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --chimeric-tag
       position: 4
   - id: exclude_chimeric
     label: exclude_chimeric
     doc: Exclude chimeric reads from the read counts
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --exclude-chimeric
       position: 4
   - id: unpaired
     label: unpaired
     doc: Allow unpaired reads to be quantified. Required for single-end libraries
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --unpaired
       position: 4
   - id: rpkm
     label: rpkm
     doc: Output gene RPKM values instead of TPMs
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --rpkm
       position: 4
   - id: coverage
     label: coverage
     doc: |-
       If this flag is provided, coverage statistics for each transcript will be written to a table. Otherwise, only summary coverage statistics are generated and added to the metrics table
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --coverage
       position: 4
   - id: coverage_mask
     label: coverage_mask
     doc: |-
       Sets how many bases at both ends of a transcript are masked out when computing per-base exon coverage. Default: 500bp
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --coverage-mask
       position: 4
   - id: detection_threshold
     label: detection_threshold
     doc: |-
       Number of counts on a gene to consider the gene 'detected'. Additionally, genes below this limit are excluded from 3' bias computation. Default: 5 reads
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --detection-threshold
       position: 4

   outputs:
   - id: out_gene_fragments
     label: out_gene_fragments
     type: File
     outputBinding:
       glob: |-
         "{output_dir}/{sample}.gene_fragments.gct".replace(/\{output_dir\}/g, inputs.output_dir).replace(/\{sample\}/g, inputs.sample)
       loadContents: false
   - id: out_gene_reads
     label: out_gene_reads
     type: File
     outputBinding:
       glob: |-
         "{output_dir}/{sample}.gene_reads.gct".replace(/\{output_dir\}/g, inputs.output_dir).replace(/\{sample\}/g, inputs.sample)
       loadContents: false
   - id: out_gene_tpm
     label: out_gene_tpm
     type: File
     outputBinding:
       glob: |-
         "{output_dir}/{sample}.gene_tpm.gct".replace(/\{output_dir\}/g, inputs.output_dir).replace(/\{sample\}/g, inputs.sample)
       loadContents: false
   - id: out_metrics_tsv
     label: out_metrics_tsv
     type: File
     outputBinding:
       glob: |-
         "{output_dir}/{sample}.metrics.tsv".replace(/\{output_dir\}/g, inputs.output_dir).replace(/\{sample\}/g, inputs.sample)
       loadContents: false
   - id: out_coverage_tsv
     label: out_coverage_tsv
     type:
     - File
     - 'null'
     outputBinding:
       glob: |-
         "{output_dir}/{sample}.coverage.tsv".replace(/\{output_dir\}/g, inputs.output_dir).replace(/\{sample\}/g, inputs.sample)
       loadContents: false
   - id: out_exon_reads
     label: out_exon_reads
     type: File
     outputBinding:
       glob: |-
         "{output_dir}/{sample}.exon_reads.gct".replace(/\{output_dir\}/g, inputs.output_dir).replace(/\{sample\}/g, inputs.sample)
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - rnaseqc
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: RNASeqQC


