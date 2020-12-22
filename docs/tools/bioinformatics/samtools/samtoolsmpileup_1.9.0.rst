:orphan:

SamTools: Mpileup
===================================

``SamToolsMpileup`` · *1 contributor · 2 versions*

Generate text pileup output for one or multiple BAM files. Each input file produces a separate group of pileup columns in the output.

Samtools mpileup can still produce VCF and BCF output (with -g or -u), but this feature is deprecated and will be removed in a future release. Please use bcftools mpileup for this instead. (Documentation on the deprecated options has been removed from this manual page, but older versions are available online at <http://www.htslib.org/doc/>.)

Note that there are two orthogonal ways to specify locations in the input file; via -r region and -l file. The former uses (and requires) an index to do random access while the latter streams through the file contents filtering out the specified regions, requiring no index. The two may be used in conjunction. For example a BED file containing locations of genes in chromosome 20 could be specified using -r 20 -l chr20.bed, meaning that the index is used to find chromosome 20 and then it is filtered for the regions listed in the bed file.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.samtools.mpileup.versions import SamToolsMpileup_1_9

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "samtoolsmpileup_step",
           SamToolsMpileup_1_9(
               bam=None,
           )
       )
       wf.output("out", source=samtoolsmpileup_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SamToolsMpileup:

.. code-block:: bash

   # user inputs
   janis inputs SamToolsMpileup > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam




5. Run SamToolsMpileup with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SamToolsMpileup





Information
------------

:ID: ``SamToolsMpileup``
:URL: `http://www.htslib.org/doc/samtools-mpileup.html <http://www.htslib.org/doc/samtools-mpileup.html>`_
:Versions: 1.9.0, 1.7.0
:Container: quay.io/biocontainers/samtools:1.9--h8571acd_11
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-05-19
:Updated: 2020-05-19


Outputs
-----------

======  ================  ===============
name    type              documentation
======  ================  ===============
out     stdout<TextFile>
======  ================  ===============


Additional configuration (inputs)
---------------------------------

======================  =================  =================  ==========  ========================================================================
name                    type               prefix               position  documentation
======================  =================  =================  ==========  ========================================================================
bam                     IndexedBam                                    10
illuminaEncoding        Optional<Boolean>  --illumina1.3+                 Assume the quality is in the Illumina 1.3+ encoding.
countOrphans            Optional<Boolean>  --count-orphans                do not discard anomalous read pairs
noBAQ                   Optional<Boolean>  --no-BAQ                       disable BAQ (per-Base Alignment Quality)
adjustMQ                Optional<Integer>  --adjust-MQ                    adjust mapping quality; recommended:50, disable:0 [0]
maxDepth                Optional<Integer>  --max-depth                    max per-file depth; avoids excessive memory usage [8000]
redoBAQ                 Optional<Boolean>  --redo-BAQ                     recalculate BAQ on the fly, ignore existing BQs
fastaRef                Optional<File>     --fasta-ref                    skip unlisted positions (chr pos) or regions (BED)
excludeRG               Optional<File>     --exclude-RG                   exclude read groups listed in FILE
positions               Optional<File>     --positions                    skip unlisted positions (chr pos) or regions (BED)
minBQ                   Optional<Integer>  --min-BQ                       Minimum base quality for a base to be considered [13]
minMQ                   Optional<Integer>  --min-MQ                       skip alignments with mapQ smaller than INT [0]
region                  Optional<String>   --region                       region in which pileup is generated
ignoreRG                Optional<Boolean>  --ignore-RG                    ignore RG tags (one BAM = one sample)
inclFlags               Optional<String>   --incl-flags                   required flags: skip reads with mask bits unset []
exclFlags               Optional<String>   --excl-flags                   filter flags: skip reads with mask bits set [UNMAP,SECONDARY,QCFAIL,DUP]
ignoreOverlaps          Optional<Boolean>  --ignore-overlaps              disable read-pair overlap detection
outputBP                Optional<Boolean>  --output-BP                    output base positions on reads
outputMQ                Optional<Boolean>  --output-MQ                    output mapping quality
outputQNAME             Optional<Boolean>  --output-QNAME                 output read names
allPositions            Optional<Boolean>  -a                             output all positions (including zero depth)
absolutelyAllPositions  Optional<Boolean>                                 output absolutely all positions, including unused ref. sequences
reference               Optional<File>     --reference                    Reference sequence FASTA FILE [null]
======================  =================  =================  ==========  ========================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task SamToolsMpileup {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Boolean? illuminaEncoding
       Boolean? countOrphans
       Boolean? noBAQ
       Int? adjustMQ
       Int? maxDepth
       Boolean? redoBAQ
       File? fastaRef
       File? excludeRG
       File? positions
       Int? minBQ
       Int? minMQ
       String? region
       Boolean? ignoreRG
       String? inclFlags
       String? exclFlags
       Boolean? ignoreOverlaps
       Boolean? outputBP
       Boolean? outputMQ
       Boolean? outputQNAME
       Boolean? allPositions
       Boolean? absolutelyAllPositions
       File? reference
       File bam
       File bam_bai
     }
     command <<<
       set -e
       samtools mpileup \
         ~{if (defined(illuminaEncoding) && select_first([illuminaEncoding])) then "--illumina1.3+" else ""} \
         ~{if (defined(countOrphans) && select_first([countOrphans])) then "--count-orphans" else ""} \
         ~{if (defined(noBAQ) && select_first([noBAQ])) then "--no-BAQ" else ""} \
         ~{if defined(adjustMQ) then ("--adjust-MQ " + adjustMQ) else ''} \
         ~{if defined(maxDepth) then ("--max-depth " + maxDepth) else ''} \
         ~{if (defined(redoBAQ) && select_first([redoBAQ])) then "--redo-BAQ" else ""} \
         ~{if defined(fastaRef) then ("--fasta-ref '" + fastaRef + "'") else ""} \
         ~{if defined(excludeRG) then ("--exclude-RG '" + excludeRG + "'") else ""} \
         ~{if defined(positions) then ("--positions '" + positions + "'") else ""} \
         ~{if defined(minBQ) then ("--min-BQ " + minBQ) else ''} \
         ~{if defined(minMQ) then ("--min-MQ " + minMQ) else ''} \
         ~{if defined(region) then ("--region '" + region + "'") else ""} \
         ~{if (defined(ignoreRG) && select_first([ignoreRG])) then "--ignore-RG" else ""} \
         ~{if defined(inclFlags) then ("--incl-flags '" + inclFlags + "'") else ""} \
         ~{if defined(exclFlags) then ("--excl-flags '" + exclFlags + "'") else ""} \
         ~{if (defined(ignoreOverlaps) && select_first([ignoreOverlaps])) then "--ignore-overlaps" else ""} \
         ~{if (defined(outputBP) && select_first([outputBP])) then "--output-BP" else ""} \
         ~{if (defined(outputMQ) && select_first([outputMQ])) then "--output-MQ" else ""} \
         ~{if (defined(outputQNAME) && select_first([outputQNAME])) then "--output-QNAME" else ""} \
         ~{if (defined(allPositions) && select_first([allPositions])) then "-a" else ""} \
         ~{if defined(reference) then ("--reference '" + reference + "'") else ""} \
         '~{bam}'
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
       File out = stdout()
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'SamTools: Mpileup'
   doc: |-
     Generate text pileup output for one or multiple BAM files. Each input file produces a separate group of pileup columns in the output.

     Samtools mpileup can still produce VCF and BCF output (with -g or -u), but this feature is deprecated and will be removed in a future release. Please use bcftools mpileup for this instead. (Documentation on the deprecated options has been removed from this manual page, but older versions are available online at <http://www.htslib.org/doc/>.)

     Note that there are two orthogonal ways to specify locations in the input file; via -r region and -l file. The former uses (and requires) an index to do random access while the latter streams through the file contents filtering out the specified regions, requiring no index. The two may be used in conjunction. For example a BED file containing locations of genes in chromosome 20 could be specified using -r 20 -l chr20.bed, meaning that the index is used to find chromosome 20 and then it is filtered for the regions listed in the bed file.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/samtools:1.9--h8571acd_11

   inputs:
   - id: illuminaEncoding
     label: illuminaEncoding
     doc: Assume the quality is in the Illumina 1.3+ encoding.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --illumina1.3+
   - id: countOrphans
     label: countOrphans
     doc: do not discard anomalous read pairs
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --count-orphans
   - id: noBAQ
     label: noBAQ
     doc: disable BAQ (per-Base Alignment Quality)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-BAQ
   - id: adjustMQ
     label: adjustMQ
     doc: adjust mapping quality; recommended:50, disable:0 [0]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --adjust-MQ
   - id: maxDepth
     label: maxDepth
     doc: max per-file depth; avoids excessive memory usage [8000]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-depth
   - id: redoBAQ
     label: redoBAQ
     doc: recalculate BAQ on the fly, ignore existing BQs
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --redo-BAQ
   - id: fastaRef
     label: fastaRef
     doc: ' skip unlisted positions (chr pos) or regions (BED)'
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --fasta-ref
   - id: excludeRG
     label: excludeRG
     doc: exclude read groups listed in FILE
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --exclude-RG
   - id: positions
     label: positions
     doc: skip unlisted positions (chr pos) or regions (BED)
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --positions
   - id: minBQ
     label: minBQ
     doc: Minimum base quality for a base to be considered [13]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-BQ
   - id: minMQ
     label: minMQ
     doc: skip alignments with mapQ smaller than INT [0]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-MQ
   - id: region
     label: region
     doc: region in which pileup is generated
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --region
   - id: ignoreRG
     label: ignoreRG
     doc: ignore RG tags (one BAM = one sample)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ignore-RG
   - id: inclFlags
     label: inclFlags
     doc: 'required flags: skip reads with mask bits unset []'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --incl-flags
   - id: exclFlags
     label: exclFlags
     doc: 'filter flags: skip reads with mask bits set [UNMAP,SECONDARY,QCFAIL,DUP]'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --excl-flags
   - id: ignoreOverlaps
     label: ignoreOverlaps
     doc: disable read-pair overlap detection
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ignore-overlaps
   - id: outputBP
     label: outputBP
     doc: output base positions on reads
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --output-BP
   - id: outputMQ
     label: outputMQ
     doc: output mapping quality
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --output-MQ
   - id: outputQNAME
     label: outputQNAME
     doc: output read names
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --output-QNAME
   - id: allPositions
     label: allPositions
     doc: output all positions (including zero depth)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -a
   - id: absolutelyAllPositions
     label: absolutelyAllPositions
     doc: output absolutely all positions, including unused ref. sequences
     type:
     - boolean
     - 'null'
   - id: reference
     label: reference
     doc: Reference sequence FASTA FILE [null]
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --reference
   - id: bam
     label: bam
     type: File
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       position: 10

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - samtools
   - mpileup
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: SamToolsMpileup


