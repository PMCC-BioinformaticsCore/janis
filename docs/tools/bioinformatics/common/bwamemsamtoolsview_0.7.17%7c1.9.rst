:orphan:

Bwa mem + Samtools View
============================================

``BwaMemSamtoolsView`` · *1 contributor · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.common.bwamem_samtoolsview import BwaMem_SamToolsView

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bwamemsamtoolsview_step",
           BwaMem_SamToolsView(
               reference=None,
               reads=None,
               sampleName=None,
           )
       )
       wf.output("out", source=bwamemsamtoolsview_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for BwaMemSamtoolsView:

.. code-block:: bash

   # user inputs
   janis inputs BwaMemSamtoolsView > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reads:
       - reads_0.fastq.gz
       - reads_1.fastq.gz
       reference: reference.fasta
       sampleName: <value>




5. Run BwaMemSamtoolsView with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       BwaMemSamtoolsView





Information
------------

:ID: ``BwaMemSamtoolsView``
:URL: *No URL to the documentation was provided*
:Versions: 0.7.17|1.9
:Container: michaelfranklin/bwasamtools:0.7.17-1.9
:Authors: Michael Franklin
:Citations: None
:Created: 2019-05-10
:Updated: 2020-07-14


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     BAM
======  ======  ===============


Additional configuration (inputs)
---------------------------------

=========================================  ========================  ============  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                                       type                      prefix          position  documentation
=========================================  ========================  ============  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
reference                                  FastaWithIndexes                                 2
reads                                      FastqGzPair                                      3
sampleName                                 String                                              Used to construct the readGroupHeaderLine with format: '@RG\tID:{name}\tSM:{name}\tLB:{name}\tPL:ILLUMINA'
mates                                      Optional<FastqGzPair>                            4
outputFilename                             Optional<Filename>        -o                     8  output file name [stdout]
platformTechnology                         Optional<String>                                    (ReadGroup: PL) Used to construct the readGroupHeaderLine, defaults: ILLUMINA
minimumSeedLength                          Optional<Integer>         -k                     2  Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. (Default: 19)
batchSize                                  Optional<Integer>         -K                     2  Process INT input bases in each batch regardless of the number of threads in use [10000000*nThreads]. By default, the batch size is proportional to the number of threads in use. Because the inferred insert size distribution slightly depends on the batch size, using different number of threads may produce different output. Specifying this option helps reproducibility.
useSoftClippingForSupplementaryAlignments  Optional<Boolean>         -Y                     2  Use soft clipping CIGAR operation for supplementary alignments. By default, BWA-MEM uses soft clipping for the primary alignment and hard clipping for supplementary alignments.
bandwidth                                  Optional<Integer>         -w                     2  Essentially, gaps longer than ${bandWidth} will not be found. Note that the maximum gap length is also affected by the scoring matrix and the hit length, not solely determined by this option. (Default: 100)
offDiagonalXDropoff                        Optional<Integer>         -d                     2  (Z-dropoff): Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. (Default: 100)
reseedTrigger                              Optional<Float>           -r                     2  Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. (Default: 1.5)
occurenceDiscard                           Optional<Integer>         -c                     2  Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. (Default: 10000)
performSW                                  Optional<Boolean>         -P                     2  In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.
matchingScore                              Optional<Integer>         -A                     2  Matching score. (Default: 1)
mismatchPenalty                            Optional<Integer>         -B                     2  Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. (Default: 4)
openGapPenalty                             Optional<Integer>         -O                     2  Gap open penalty. (Default: 6)
gapExtensionPenalty                        Optional<Integer>         -E                     2  Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). (Default: 1)
clippingPenalty                            Optional<Integer>         -L                     2  Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. (Default: 5)
unpairedReadPenalty                        Optional<Integer>         -U                     2  Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. It compares these two scores to determine whether we should force pairing. (Default: 9)
assumeInterleavedFirstInput                Optional<Boolean>         -p                     2  Assume the first input query file is interleaved paired-end FASTA/Q.
outputAlignmentThreshold                   Optional<Integer>         -T                     2  Don’t output alignment with score lower than INT. Only affects output. (Default: 30)
outputAllElements                          Optional<Boolean>         -a                     2  Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.
appendComments                             Optional<Boolean>         -C                     2  Append append FASTA/Q comment to SAM output. This option can be used to transfer read meta information (e.g. barcode) to the SAM output. Note that the FASTA/Q comment (the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). Malformated comments lead to incorrect SAM output.
hardClipping                               Optional<Boolean>         -H                     2  Use hard clipping ’H’ in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences.
markShorterSplits                          Optional<Boolean>         -M                     2  Mark shorter split hits as secondary (for Picard compatibility).
verboseLevel                               Optional<Integer>         -v                     2  Control the verbose level of the output. This option has not been fully supported throughout BWA. Ideally, a value: 0 for disabling all the output to stderr; 1 for outputting errors only; 2 for warnings and errors; 3 for all normal messages; 4 or higher for debugging. When this option takes value 4, the output is not SAM. (Default: 3)
skippedReadsOutputFilename                 Optional<String>          -U                     8  output reads not selected by filters to FILE [null]
referenceIndex                             Optional<File>            -t                     8  FILE listing reference names and lengths (see long help) [null]
intervals                                  Optional<bed>             -L                     8  only include reads overlapping this BED FILE [null]
includeReadsInReadGroup                    Optional<String>          -r                     8  only include reads in read group STR [null]
includeReadsInFile                         Optional<File>            -R                     8  only include reads with read group listed in FILE [null]
includeReadsWithQuality                    Optional<Integer>         -q                     8  only include reads with mapping quality >= INT [0]
includeReadsInLibrary                      Optional<String>          -l                     8  only include reads in library STR [null]
includeReadsWithCIGAROps                   Optional<Integer>         -m                     8  only include reads with number of CIGAR operations consuming query sequence >= INT [0]
includeReadsWithAllFLAGs                   Optional<Array<Integer>>  -f                     8  only include reads with all of the FLAGs in INT present [0]
includeReadsWithoutFLAGs                   Optional<Array<Integer>>  -F                     8  only include reads with none of the FLAGS in INT present [0]
excludeReadsWithAllFLAGs                   Optional<Array<Integer>>  -G                     8  only EXCLUDE reads with all of the FLAGs in INT present [0] fraction of templates/read pairs to keep; INT part sets seed)
useMultiRegionIterator                     Optional<Boolean>         -M                     8  use the multi-region iterator (increases the speed, removes duplicates and outputs the reads as they are ordered in the file)
readTagToStrip                             Optional<String>          -x                     8  read tag to strip (repeatable) [null]
collapseBackwardCIGAROps                   Optional<Boolean>         -B                     8  collapse the backward CIGAR operation Specify a single input file format option in the form of OPTION or OPTION=VALUE
outputFmt                                  Optional<String>          --output-fmt           8  (OPT[, -O)  Specify output format (SAM, BAM, CRAM) Specify a single output file format option in the form of OPTION or OPTION=VALUE
=========================================  ========================  ============  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task BwaMemSamtoolsView {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       Array[File] reads
       Array[File]? mates
       String? outputFilename
       String sampleName
       String? platformTechnology
       Int? minimumSeedLength
       Int? batchSize
       Boolean? useSoftClippingForSupplementaryAlignments
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
       Int? outputAlignmentThreshold
       Boolean? outputAllElements
       Boolean? appendComments
       Boolean? hardClipping
       Boolean? markShorterSplits
       Int? verboseLevel
       String? skippedReadsOutputFilename
       File? referenceIndex
       File? intervals
       String? includeReadsInReadGroup
       File? includeReadsInFile
       Int? includeReadsWithQuality
       String? includeReadsInLibrary
       Int? includeReadsWithCIGAROps
       Array[Int]? includeReadsWithAllFLAGs
       Array[Int]? includeReadsWithoutFLAGs
       Array[Int]? excludeReadsWithAllFLAGs
       Boolean? useMultiRegionIterator
       String? readTagToStrip
       Boolean? collapseBackwardCIGAROps
       String? outputFmt
     }
     command <<<
       set -e
        \
         bwa \
         mem \
         ~{reference} \
         ~{if defined(minimumSeedLength) then ("-k " + minimumSeedLength) else ''} \
         ~{if defined(batchSize) then ("-K " + batchSize) else ''} \
         ~{if (defined(useSoftClippingForSupplementaryAlignments) && select_first([useSoftClippingForSupplementaryAlignments])) then "-Y" else ""} \
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
         ~{if defined(outputAlignmentThreshold) then ("-T " + outputAlignmentThreshold) else ''} \
         ~{if (defined(outputAllElements) && select_first([outputAllElements])) then "-a" else ""} \
         ~{if (defined(appendComments) && select_first([appendComments])) then "-C" else ""} \
         ~{if (defined(hardClipping) && select_first([hardClipping])) then "-H" else ""} \
         ~{if (defined(markShorterSplits) && select_first([markShorterSplits])) then "-M" else ""} \
         ~{if defined(verboseLevel) then ("-v " + verboseLevel) else ''} \
         -R '@RG\tID:~{sampleName}\tSM:~{sampleName}\tLB:~{sampleName}\tPL:~{select_first([platformTechnology, "ILLUMINA"])}' \
         -t ~{select_first([runtime_cpu, 16, 1])} \
         ~{if length(reads) > 0 then sep(" ", reads) else ""} \
         ~{if (defined(mates) && length(select_first([mates])) > 0) then sep(" ", select_first([mates])) else ""} \
         | \
         samtools \
         view \
         -o ~{select_first([outputFilename, "~{sampleName}.bam"])} \
         ~{if defined(skippedReadsOutputFilename) then ("-U " + skippedReadsOutputFilename) else ''} \
         ~{if defined(referenceIndex) then ("-t " + referenceIndex) else ''} \
         ~{if defined(intervals) then ("-L " + intervals) else ''} \
         ~{if defined(includeReadsInReadGroup) then ("-r " + includeReadsInReadGroup) else ''} \
         ~{if defined(includeReadsInFile) then ("-R " + includeReadsInFile) else ''} \
         ~{if defined(includeReadsWithQuality) then ("-q " + includeReadsWithQuality) else ''} \
         ~{if defined(includeReadsInLibrary) then ("-l " + includeReadsInLibrary) else ''} \
         ~{if defined(includeReadsWithCIGAROps) then ("-m " + includeReadsWithCIGAROps) else ''} \
         ~{if (defined(includeReadsWithAllFLAGs) && length(select_first([includeReadsWithAllFLAGs])) > 0) then "-f " + sep(" ", select_first([includeReadsWithAllFLAGs])) else ""} \
         ~{if (defined(includeReadsWithoutFLAGs) && length(select_first([includeReadsWithoutFLAGs])) > 0) then "-F " + sep(" ", select_first([includeReadsWithoutFLAGs])) else ""} \
         ~{if (defined(excludeReadsWithAllFLAGs) && length(select_first([excludeReadsWithAllFLAGs])) > 0) then "-G " + sep(" ", select_first([excludeReadsWithAllFLAGs])) else ""} \
         ~{if (defined(useMultiRegionIterator) && select_first([useMultiRegionIterator])) then "-M" else ""} \
         ~{if defined(readTagToStrip) then ("-x " + readTagToStrip) else ''} \
         ~{if (defined(collapseBackwardCIGAROps) && select_first([collapseBackwardCIGAROps])) then "-B" else ""} \
         ~{if defined(outputFmt) then ("--output-fmt " + outputFmt) else ''} \
         -T ~{reference} \
         --threads ~{select_first([runtime_cpu, 16, 1])} \
         -h \
         -b
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 16, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "michaelfranklin/bwasamtools:0.7.17-1.9"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 16, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "~{sampleName}.bam"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Bwa mem + Samtools View

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: michaelfranklin/bwasamtools:0.7.17-1.9

   inputs:
   - id: reference
     label: reference
     type: File
     secondaryFiles:
     - pattern: .fai
     - pattern: .amb
     - pattern: .ann
     - pattern: .bwt
     - pattern: .pac
     - pattern: .sa
     - pattern: ^.dict
     inputBinding:
       position: 2
       shellQuote: false
   - id: reads
     label: reads
     type:
       type: array
       items: File
     inputBinding:
       position: 3
       shellQuote: false
   - id: mates
     label: mates
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       position: 4
       itemSeparator: ' '
       shellQuote: false
   - id: outputFilename
     label: outputFilename
     doc: output file name [stdout]
     type:
     - string
     - 'null'
     default: generated.bam
     inputBinding:
       prefix: -o
       position: 8
       valueFrom: $(inputs.sampleName).bam
       shellQuote: false
   - id: sampleName
     label: sampleName
     doc: |-
       Used to construct the readGroupHeaderLine with format: '@RG\tID:{name}\tSM:{name}\tLB:{name}\tPL:ILLUMINA'
     type: string
   - id: platformTechnology
     label: platformTechnology
     doc: '(ReadGroup: PL) Used to construct the readGroupHeaderLine, defaults: ILLUMINA'
     type: string
     default: ILLUMINA
   - id: minimumSeedLength
     label: minimumSeedLength
     doc: |-
       Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. (Default: 19)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -k
       position: 2
       shellQuote: false
   - id: batchSize
     label: batchSize
     doc: |-
       Process INT input bases in each batch regardless of the number of threads in use [10000000*nThreads]. By default, the batch size is proportional to the number of threads in use. Because the inferred insert size distribution slightly depends on the batch size, using different number of threads may produce different output. Specifying this option helps reproducibility.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -K
       position: 2
       shellQuote: false
   - id: useSoftClippingForSupplementaryAlignments
     label: useSoftClippingForSupplementaryAlignments
     doc: |-
       Use soft clipping CIGAR operation for supplementary alignments. By default, BWA-MEM uses soft clipping for the primary alignment and hard clipping for supplementary alignments.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -Y
       position: 2
       shellQuote: false
   - id: bandwidth
     label: bandwidth
     doc: |-
       Essentially, gaps longer than ${bandWidth} will not be found. Note that the maximum gap length is also affected by the scoring matrix and the hit length, not solely determined by this option. (Default: 100)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -w
       position: 2
       shellQuote: false
   - id: offDiagonalXDropoff
     label: offDiagonalXDropoff
     doc: |-
       (Z-dropoff): Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. (Default: 100)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -d
       position: 2
       shellQuote: false
   - id: reseedTrigger
     label: reseedTrigger
     doc: |-
       Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. (Default: 1.5)
     type:
     - float
     - 'null'
     inputBinding:
       prefix: -r
       position: 2
       shellQuote: false
   - id: occurenceDiscard
     label: occurenceDiscard
     doc: |-
       Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. (Default: 10000)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -c
       position: 2
       shellQuote: false
   - id: performSW
     label: performSW
     doc: |-
       In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -P
       position: 2
       shellQuote: false
   - id: matchingScore
     label: matchingScore
     doc: 'Matching score. (Default: 1)'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -A
       position: 2
       shellQuote: false
   - id: mismatchPenalty
     label: mismatchPenalty
     doc: |-
       Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. (Default: 4)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -B
       position: 2
       shellQuote: false
   - id: openGapPenalty
     label: openGapPenalty
     doc: 'Gap open penalty. (Default: 6)'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -O
       position: 2
       shellQuote: false
   - id: gapExtensionPenalty
     label: gapExtensionPenalty
     doc: |-
       Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). (Default: 1)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -E
       position: 2
       shellQuote: false
   - id: clippingPenalty
     label: clippingPenalty
     doc: |-
       Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. (Default: 5)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -L
       position: 2
       shellQuote: false
   - id: unpairedReadPenalty
     label: unpairedReadPenalty
     doc: |-
       Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. It compares these two scores to determine whether we should force pairing. (Default: 9)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -U
       position: 2
       shellQuote: false
   - id: assumeInterleavedFirstInput
     label: assumeInterleavedFirstInput
     doc: 'Assume the first input query file is interleaved paired-end FASTA/Q. '
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -p
       position: 2
       shellQuote: false
   - id: outputAlignmentThreshold
     label: outputAlignmentThreshold
     doc: |-
       Don’t output alignment with score lower than INT. Only affects output. (Default: 30)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -T
       position: 2
       shellQuote: false
   - id: outputAllElements
     label: outputAllElements
     doc: |-
       Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -a
       position: 2
       shellQuote: false
   - id: appendComments
     label: appendComments
     doc: |-
       Append append FASTA/Q comment to SAM output. This option can be used to transfer read meta information (e.g. barcode) to the SAM output. Note that the FASTA/Q comment (the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). Malformated comments lead to incorrect SAM output.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -C
       position: 2
       shellQuote: false
   - id: hardClipping
     label: hardClipping
     doc: |-
       Use hard clipping ’H’ in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -H
       position: 2
       shellQuote: false
   - id: markShorterSplits
     label: markShorterSplits
     doc: Mark shorter split hits as secondary (for Picard compatibility).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -M
       position: 2
       shellQuote: false
   - id: verboseLevel
     label: verboseLevel
     doc: |-
       Control the verbose level of the output. This option has not been fully supported throughout BWA. Ideally, a value: 0 for disabling all the output to stderr; 1 for outputting errors only; 2 for warnings and errors; 3 for all normal messages; 4 or higher for debugging. When this option takes value 4, the output is not SAM. (Default: 3)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -v
       position: 2
       shellQuote: false
   - id: skippedReadsOutputFilename
     label: skippedReadsOutputFilename
     doc: output reads not selected by filters to FILE [null]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -U
       position: 8
       shellQuote: false
   - id: referenceIndex
     label: referenceIndex
     doc: FILE listing reference names and lengths (see long help) [null]
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -t
       position: 8
       shellQuote: false
   - id: intervals
     label: intervals
     doc: only include reads overlapping this BED FILE [null]
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -L
       position: 8
       shellQuote: false
   - id: includeReadsInReadGroup
     label: includeReadsInReadGroup
     doc: only include reads in read group STR [null]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -r
       position: 8
       shellQuote: false
   - id: includeReadsInFile
     label: includeReadsInFile
     doc: only include reads with read group listed in FILE [null]
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -R
       position: 8
       shellQuote: false
   - id: includeReadsWithQuality
     label: includeReadsWithQuality
     doc: only include reads with mapping quality >= INT [0]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -q
       position: 8
       shellQuote: false
   - id: includeReadsInLibrary
     label: includeReadsInLibrary
     doc: only include reads in library STR [null]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -l
       position: 8
       shellQuote: false
   - id: includeReadsWithCIGAROps
     label: includeReadsWithCIGAROps
     doc: |-
       only include reads with number of CIGAR operations consuming query sequence >= INT [0]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -m
       position: 8
       shellQuote: false
   - id: includeReadsWithAllFLAGs
     label: includeReadsWithAllFLAGs
     doc: only include reads with all of the FLAGs in INT present [0]
     type:
     - type: array
       items: int
     - 'null'
     inputBinding:
       prefix: -f
       position: 8
       itemSeparator: ' '
       shellQuote: false
   - id: includeReadsWithoutFLAGs
     label: includeReadsWithoutFLAGs
     doc: only include reads with none of the FLAGS in INT present [0]
     type:
     - type: array
       items: int
     - 'null'
     inputBinding:
       prefix: -F
       position: 8
       itemSeparator: ' '
       shellQuote: false
   - id: excludeReadsWithAllFLAGs
     label: excludeReadsWithAllFLAGs
     doc: |-
       only EXCLUDE reads with all of the FLAGs in INT present [0] fraction of templates/read pairs to keep; INT part sets seed)
     type:
     - type: array
       items: int
     - 'null'
     inputBinding:
       prefix: -G
       position: 8
       itemSeparator: ' '
       shellQuote: false
   - id: useMultiRegionIterator
     label: useMultiRegionIterator
     doc: |-
       use the multi-region iterator (increases the speed, removes duplicates and outputs the reads as they are ordered in the file)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -M
       position: 8
       shellQuote: false
   - id: readTagToStrip
     label: readTagToStrip
     doc: read tag to strip (repeatable) [null]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -x
       position: 8
       shellQuote: false
   - id: collapseBackwardCIGAROps
     label: collapseBackwardCIGAROps
     doc: |-
       collapse the backward CIGAR operation Specify a single input file format option in the form of OPTION or OPTION=VALUE
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -B
       position: 8
       shellQuote: false
   - id: outputFmt
     label: outputFmt
     doc: |-
       (OPT[, -O)  Specify output format (SAM, BAM, CRAM) Specify a single output file format option in the form of OPTION or OPTION=VALUE
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --output-fmt
       position: 8
       shellQuote: false

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: $(inputs.sampleName).bam
       loadContents: false
   stdout: _stdout
   stderr: _stderr
   arguments:
   - position: 0
     valueFrom: bwa
     shellQuote: false
   - position: 1
     valueFrom: mem
     shellQuote: false
   - position: 5
     valueFrom: '|'
     shellQuote: false
   - position: 6
     valueFrom: samtools
     shellQuote: false
   - position: 7
     valueFrom: view
     shellQuote: false
   - prefix: -T
     position: 8
     valueFrom: $(inputs.reference)
     shellQuote: false
   - prefix: --threads
     position: 8
     valueFrom: |-
       $([inputs.runtime_cpu, 16, 1].filter(function (inner) { return inner != null })[0])
     shellQuote: false
   - position: 8
     valueFrom: -h
     shellQuote: false
   - position: 8
     valueFrom: -b
     shellQuote: false
   - prefix: -R
     position: 2
     valueFrom: |-
       $("@RG\\tID:{name}\\tSM:{name}\\tLB:{name}\\tPL:{pl}".replace(/\{name\}/g, inputs.sampleName).replace(/\{pl\}/g, inputs.platformTechnology))
   - prefix: -t
     position: 2
     valueFrom: |-
       $([inputs.runtime_cpu, 16, 1].filter(function (inner) { return inner != null })[0])
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: BwaMemSamtoolsView


