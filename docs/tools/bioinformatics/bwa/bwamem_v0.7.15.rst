:orphan:

BWA-MEM
================

``bwamem`` · *1 contributor · 1 version*

bwa - Burrows-Wheeler Alignment Tool
BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human 
genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for 
Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. 
BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is 
the latest, is generally recommended for high-quality queries as it is faster and more accurate. 
BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

Align 70bp-1Mbp query sequences with the BWA-MEM algorithm. Briefly, the algorithm works by seeding alignments 
with maximal exact matches (MEMs) and then extending seeds with the affine-gap Smith-Waterman algorithm (SW).

If mates.fq file is absent and option -p is not set, this command regards input reads are single-end. If 'mates.fq' 
is present, this command assumes the i-th read in reads.fq and the i-th read in mates.fq constitute a read pair. 
If -p is used, the command assumes the 2i-th and the (2i+1)-th read in reads.fq constitute a read pair (such input 
file is said to be interleaved). In this case, mates.fq is ignored. In the paired-end mode, the mem command will 
infer the read orientation and the insert size distribution from a batch of reads.

The BWA-MEM algorithm performs local alignment. It may produce multiple primary alignments for different part of a 
query sequence. This is a crucial feature for long sequences. However, some tools such as Picard’s markDuplicates 
does not work with split alignments. One may consider to use option -M to flag shorter split hits as secondary.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.bwa.mem.versions import BwaMem_0_7_15

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bwamem_step",
           BwaMem_0_7_15(
               reference=None,
               reads=None,
           )
       )
       wf.output("out", source=bwamem_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bwamem:

.. code-block:: bash

   # user inputs
   janis inputs bwamem > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reads:
       - reads_0.fastq.gz
       - reads_1.fastq.gz
       reference: reference.fasta




5. Run bwamem with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bwamem





Information
------------

:ID: ``bwamem``
:URL: `http://bio-bwa.sourceforge.net/bwa.shtml#3 <http://bio-bwa.sourceforge.net/bwa.shtml#3>`_
:Versions: v0.7.15
:Container: biocontainers/bwa:v0.7.15_cv3
:Authors: Michael Franklin
:Citations: The BWA-MEM algorithm has not been published yet.
:Created: 2018-12-24
:Updated: 2019-07-23


Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     stdout<SAM>
======  ===========  ===============


Additional configuration (inputs)
---------------------------------

===========================  =====================  ========  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                         type                   prefix      position  documentation
===========================  =====================  ========  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
reference                    FastaBwa                                  9
reads                        Array<FastqGz>                           10
mates                        Optional<FastqGzPair>                    11
outputFilename               Optional<Filename>
threads                      Optional<Integer>      -t                    Number of threads. (default = 1)
minimumSeedLength            Optional<Integer>      -k                    Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. (Default: 19)
bandwidth                    Optional<Integer>      -w                    Essentially, gaps longer than ${bandWidth} will not be found. Note that the maximum gap length is also affected by the scoring matrix and the hit length, not solely determined by this option. (Default: 100)
offDiagonalXDropoff          Optional<Integer>      -d                    (Z-dropoff): Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. (Default: 100)
reseedTrigger                Optional<Float>        -r                    Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. (Default: 1.5)
occurenceDiscard             Optional<Integer>      -c                    Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. (Default: 10000)
performSW                    Optional<Boolean>      -P                    In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.
matchingScore                Optional<Integer>      -A                    Matching score. (Default: 1)
mismatchPenalty              Optional<Integer>      -B                    Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. (Default: 4)
openGapPenalty               Optional<Integer>      -O                    Gap open penalty. (Default: 6)
gapExtensionPenalty          Optional<Integer>      -E                    Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). (Default: 1)
clippingPenalty              Optional<Integer>      -L                    Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. (Default: 5)
unpairedReadPenalty          Optional<Integer>      -U                    Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. It compares these two scores to determine whether we should force pairing. (Default: 9)
assumeInterleavedFirstInput  Optional<Boolean>      -p                    Assume the first input query file is interleaved paired-end FASTA/Q.
readGroupHeaderLine          Optional<String>       -R                    Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’. (Default=null)
outputAlignmentThreshold     Optional<Integer>      -T                    Don’t output alignment with score lower than INT. Only affects output. (Default: 30)
outputAllElements            Optional<Boolean>      -a                    Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.
appendComments               Optional<Boolean>      -C                    Append append FASTA/Q comment to SAM output. This option can be used to transfer read meta information (e.g. barcode) to the SAM output. Note that the FASTA/Q comment (the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). Malformated comments lead to incorrect SAM output.
hardClipping                 Optional<Boolean>      -H                    Use hard clipping ’H’ in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences.
markShorterSplits            Optional<Boolean>      -M                    Mark shorter split hits as secondary (for Picard compatibility).
verboseLevel                 Optional<Integer>      -v                    Control the verbose level of the output. This option has not been fully supported throughout BWA. Ideally, a value: 0 for disabling all the output to stderr; 1 for outputting errors only; 2 for warnings and errors; 3 for all normal messages; 4 or higher for debugging. When this option takes value 4, the output is not SAM. (Default: 3)
===========================  =====================  ========  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task bwamem {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File reference
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       Array[File] reads
       Array[File]? mates
       String? outputFilename
       Int? threads
       Int? minimumSeedLength
       Int? bandwidth
       Int? offDiagonalXDropoff
       Float? reseedTrigger
       Int? occurenceDiscard
       Boolean? performSW
       Int? matchingScore
       Int? mismatchPenalty
       Int? openGapPenalty
       Int? gapExtensionPenalty
       Int? clippingPenalty
       Int? unpairedReadPenalty
       Boolean? assumeInterleavedFirstInput
       String? readGroupHeaderLine
       Int? outputAlignmentThreshold
       Boolean? outputAllElements
       Boolean? appendComments
       Boolean? hardClipping
       Boolean? markShorterSplits
       Int? verboseLevel
     }
     command <<<
       set -e
       bwa mem \
         ~{if defined(select_first([threads, select_first([runtime_cpu, 1])])) then ("-t " + select_first([threads, select_first([runtime_cpu, 1])])) else ''} \
         ~{if defined(minimumSeedLength) then ("-k " + minimumSeedLength) else ''} \
         ~{if defined(bandwidth) then ("-w " + bandwidth) else ''} \
         ~{if defined(offDiagonalXDropoff) then ("-d " + offDiagonalXDropoff) else ''} \
         ~{if defined(reseedTrigger) then ("-r " + reseedTrigger) else ''} \
         ~{if defined(occurenceDiscard) then ("-c " + occurenceDiscard) else ''} \
         ~{if (defined(performSW) && select_first([performSW])) then "-P" else ""} \
         ~{if defined(matchingScore) then ("-A " + matchingScore) else ''} \
         ~{if defined(mismatchPenalty) then ("-B " + mismatchPenalty) else ''} \
         ~{if defined(openGapPenalty) then ("-O " + openGapPenalty) else ''} \
         ~{if defined(gapExtensionPenalty) then ("-E " + gapExtensionPenalty) else ''} \
         ~{if defined(clippingPenalty) then ("-L " + clippingPenalty) else ''} \
         ~{if defined(unpairedReadPenalty) then ("-U " + unpairedReadPenalty) else ''} \
         ~{if (defined(assumeInterleavedFirstInput) && select_first([assumeInterleavedFirstInput])) then "-p" else ""} \
         ~{if defined(readGroupHeaderLine) then ("-R '" + readGroupHeaderLine + "'") else ""} \
         ~{if defined(outputAlignmentThreshold) then ("-T " + outputAlignmentThreshold) else ''} \
         ~{if (defined(outputAllElements) && select_first([outputAllElements])) then "-a" else ""} \
         ~{if (defined(appendComments) && select_first([appendComments])) then "-C" else ""} \
         ~{if (defined(hardClipping) && select_first([hardClipping])) then "-H" else ""} \
         ~{if (defined(markShorterSplits) && select_first([markShorterSplits])) then "-M" else ""} \
         ~{if defined(verboseLevel) then ("-v " + verboseLevel) else ''} \
         '~{reference}' \
         ~{if length(reads) > 0 then "'" + sep("' '", reads) + "'" else ""} \
         ~{if (defined(mates) && length(select_first([mates])) > 0) then "'" + sep("' '", select_first([mates])) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 16, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "biocontainers/bwa:v0.7.15_cv3"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 16, 4])}G"
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
   label: BWA-MEM
   doc: |-
     bwa - Burrows-Wheeler Alignment Tool
     BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human 
     genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for 
     Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. 
     BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is 
     the latest, is generally recommended for high-quality queries as it is faster and more accurate. 
     BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

     Align 70bp-1Mbp query sequences with the BWA-MEM algorithm. Briefly, the algorithm works by seeding alignments 
     with maximal exact matches (MEMs) and then extending seeds with the affine-gap Smith-Waterman algorithm (SW).

     If mates.fq file is absent and option -p is not set, this command regards input reads are single-end. If 'mates.fq' 
     is present, this command assumes the i-th read in reads.fq and the i-th read in mates.fq constitute a read pair. 
     If -p is used, the command assumes the 2i-th and the (2i+1)-th read in reads.fq constitute a read pair (such input 
     file is said to be interleaved). In this case, mates.fq is ignored. In the paired-end mode, the mem command will 
     infer the read orientation and the insert size distribution from a batch of reads.

     The BWA-MEM algorithm performs local alignment. It may produce multiple primary alignments for different part of a 
     query sequence. This is a crucial feature for long sequences. However, some tools such as Picard’s markDuplicates 
     does not work with split alignments. One may consider to use option -M to flag shorter split hits as secondary.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: biocontainers/bwa:v0.7.15_cv3

   inputs:
   - id: reference
     label: reference
     type: File
     secondaryFiles:
     - pattern: .amb
     - pattern: .ann
     - pattern: .bwt
     - pattern: .pac
     - pattern: .sa
     inputBinding:
       position: 9
   - id: reads
     label: reads
     type:
       type: array
       items: File
     inputBinding:
       position: 10
   - id: mates
     label: mates
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       position: 11
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.sam
   - id: threads
     label: threads
     doc: Number of threads. (default = 1)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -t
       valueFrom: |-
         $([inputs.runtime_cpu, 16, 1].filter(function (inner) { return inner != null })[0])
   - id: minimumSeedLength
     label: minimumSeedLength
     doc: |-
       Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. (Default: 19)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -k
   - id: bandwidth
     label: bandwidth
     doc: |-
       Essentially, gaps longer than ${bandWidth} will not be found. Note that the maximum gap length is also affected by the scoring matrix and the hit length, not solely determined by this option. (Default: 100)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -w
   - id: offDiagonalXDropoff
     label: offDiagonalXDropoff
     doc: |-
       (Z-dropoff): Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. (Default: 100)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -d
   - id: reseedTrigger
     label: reseedTrigger
     doc: |-
       Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. (Default: 1.5)
     type:
     - float
     - 'null'
     inputBinding:
       prefix: -r
   - id: occurenceDiscard
     label: occurenceDiscard
     doc: |-
       Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. (Default: 10000)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -c
   - id: performSW
     label: performSW
     doc: |-
       In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -P
   - id: matchingScore
     label: matchingScore
     doc: 'Matching score. (Default: 1)'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -A
   - id: mismatchPenalty
     label: mismatchPenalty
     doc: |-
       Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. (Default: 4)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -B
   - id: openGapPenalty
     label: openGapPenalty
     doc: 'Gap open penalty. (Default: 6)'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -O
   - id: gapExtensionPenalty
     label: gapExtensionPenalty
     doc: |-
       Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). (Default: 1)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -E
   - id: clippingPenalty
     label: clippingPenalty
     doc: |-
       Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. (Default: 5)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -L
   - id: unpairedReadPenalty
     label: unpairedReadPenalty
     doc: |-
       Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. It compares these two scores to determine whether we should force pairing. (Default: 9)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -U
   - id: assumeInterleavedFirstInput
     label: assumeInterleavedFirstInput
     doc: 'Assume the first input query file is interleaved paired-end FASTA/Q. '
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -p
   - id: readGroupHeaderLine
     label: readGroupHeaderLine
     doc: |-
       Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’. (Default=null)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -R
   - id: outputAlignmentThreshold
     label: outputAlignmentThreshold
     doc: |-
       Don’t output alignment with score lower than INT. Only affects output. (Default: 30)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -T
   - id: outputAllElements
     label: outputAllElements
     doc: |-
       Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -a
   - id: appendComments
     label: appendComments
     doc: |-
       Append append FASTA/Q comment to SAM output. This option can be used to transfer read meta information (e.g. barcode) to the SAM output. Note that the FASTA/Q comment (the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). Malformated comments lead to incorrect SAM output.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -C
   - id: hardClipping
     label: hardClipping
     doc: |-
       Use hard clipping ’H’ in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -H
   - id: markShorterSplits
     label: markShorterSplits
     doc: Mark shorter split hits as secondary (for Picard compatibility).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -M
   - id: verboseLevel
     label: verboseLevel
     doc: |-
       Control the verbose level of the output. This option has not been fully supported throughout BWA. Ideally, a value: 0 for disabling all the output to stderr; 1 for outputting errors only; 2 for warnings and errors; 3 for all normal messages; 4 or higher for debugging. When this option takes value 4, the output is not SAM. (Default: 3)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -v

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - bwa
   - mem
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: bwamem


