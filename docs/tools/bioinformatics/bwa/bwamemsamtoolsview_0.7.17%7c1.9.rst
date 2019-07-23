:orphan:


Bwa mem + Samtools View
============================================
Tool identifier: ``BwaMemSamtoolsView``
Tool path: ``janis_bioinformatics.tools.common.bwamem_samtoolsview import BwaMem_SamToolsView``

Version: 0.7.17|1.9
Docker: ``michaelfranklin/bwasamtools:0.7.17-1.9``



Documentation
-------------

URL
******
*No URL to the documentation was provided*

Description
*********
*No documentation was provided: `contribute one <https://github.com/illusional>`_*

Outputs
-------
======  ======  ===============
name    type    documentation
======  ======  ===============
out     BAM
======  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

==========  =============  ========  ==========  ==========================================================================================================
name        type           prefix      position  documentation
==========  =============  ========  ==========  ==========================================================================================================
reference   FastaWithDict                     2
reads       Fastq                             3
sampleName  String                               Used to construct the readGroupHeaderLine with format: '@RG\tID:{name}\tSM:{name}\tLB:{name}\tPL:ILLUMINA'
==========  =============  ========  ==========  ==========================================================================================================

Optional inputs
***************

===========================  ========================  ============  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                         type                      prefix          position  documentation
===========================  ========================  ============  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
mates                        Optional<Fastq>                                  4
outputFilename               Optional<Filename>        -o                     8  output file name [stdout]
minimumSeedLength            Optional<Integer>         -k                     2  Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. (Default: 19)
bandwidth                    Optional<Integer>         -w                     2  Essentially, gaps longer than ${bandWidth} will not be found. Note that the maximum gap length is also affected by the scoring matrix and the hit length, not solely determined by this option. (Default: 100)
offDiagonalXDropoff          Optional<Integer>         -d                     2  (Z-dropoff): Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. (Default: 100)
reseedTrigger                Optional<Float>           -r                     2  Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. (Default: 1.5)
occurenceDiscard             Optional<Integer>         -c                     2  Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. (Default: 10000)
performSW                    Optional<Boolean>         -P                     2  In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.
matchingScore                Optional<Integer>         -A                     2  Matching score. (Default: 1)
mismatchPenalty              Optional<Integer>         -B                     2  Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. (Default: 4)
openGapPenalty               Optional<Integer>         -O                     2  Gap open penalty. (Default: 6)
gapExtensionPenalty          Optional<Integer>         -E                     2  Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). (Default: 1)
clippingPenalty              Optional<Integer>         -L                     2  Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. (Default: 5)
unpairedReadPenalty          Optional<Integer>         -U                     2  Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. It compares these two scores to determine whether we should force pairing. (Default: 9)
assumeInterleavedFirstInput  Optional<Boolean>         -p                     2  Assume the first input query file is interleaved paired-end FASTA/Q.
outputAlignmentThreshold     Optional<Integer>         -T                     2  Don’t output alignment with score lower than INT. Only affects output. (Default: 30)
outputAllElements            Optional<Boolean>         -a                     2  Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.
appendComments               Optional<Boolean>         -C                     2  Append append FASTA/Q comment to SAM output. This option can be used to transfer read meta information (e.g. barcode) to the SAM output. Note that the FASTA/Q comment (the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). Malformated comments lead to incorrect SAM output.
hardClipping                 Optional<Boolean>         -H                     2  Use hard clipping ’H’ in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences.
markShorterSplits            Optional<Boolean>         -M                     2  Mark shorter split hits as secondary (for Picard compatibility).
verboseLevel                 Optional<Integer>         -v                     2  Control the verbose level of the output. This option has not been fully supported throughout BWA. Ideally, a value: 0 for disabling all the output to stderr; 1 for outputting errors only; 2 for warnings and errors; 3 for all normal messages; 4 or higher for debugging. When this option takes value 4, the output is not SAM. (Default: 3)
skippedReadsOutputFilename   Optional<String>          -U                     8  output reads not selected by filters to FILE [null]
referenceIndex               Optional<File>            -t                     8  FILE listing reference names and lengths (see long help) [null]
intervals                    Optional<bed>             -L                     8  only include reads overlapping this BED FILE [null]
includeReadsInReadGroup      Optional<String>          -r                     8  only include reads in read group STR [null]
includeReadsInFile           Optional<File>            -R                     8  only include reads with read group listed in FILE [null]
includeReadsWithQuality      Optional<Integer>         -q                     8  only include reads with mapping quality >= INT [0]
includeReadsInLibrary        Optional<String>          -l                     8  only include reads in library STR [null]
includeReadsWithCIGAROps     Optional<Integer>         -m                     8  only include reads with number of CIGAR operations consuming query sequence >= INT [0]
includeReadsWithAllFLAGs     Optional<Array<Integer>>  -f                     8  only include reads with all of the FLAGs in INT present [0]
includeReadsWithoutFLAGs     Optional<Array<Integer>>  -F                     8  only include reads with none of the FLAGS in INT present [0]
excludeReadsWithAllFLAGs     Optional<Array<Integer>>  -G                     8  only EXCLUDE reads with all of the FLAGs in INT present [0] fraction of templates/read pairs to keep; INT part sets seed)
useMultiRegionIterator       Optional<Boolean>         -M                     8  use the multi-region iterator (increases the speed, removes duplicates and outputs the reads as they are ordered in the file)
readTagToStrip               Optional<String>          -x                     8  read tag to strip (repeatable) [null]
collapseBackwardCIGAROps     Optional<Boolean>         -B                     8  collapse the backward CIGAR operation Specify a single input file format option in the form of OPTION or OPTION=VALUE
outputFmt                    Optional<String>          --output-fmt           8  (OPT[, -O)  Specify output format (SAM, BAM, CRAM) Specify a single output file format option in the form of OPTION or OPTION=VALUE
===========================  ========================  ============  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: **Unknown**


*Bwa mem + Samtools View was last updated on **Unknown***.
*This page was automatically generated on 2019-07-23*.
