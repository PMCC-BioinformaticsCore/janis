
BWA-MEM
================
Tool identifier: ``bwamem``

Tool path: ``from janis_bioinformatics.tools.bwa import BwaMem_0_7_15``

Documentation
-------------

Docker
******
``biocontainers/bwa:v0.7.15_cv3``

URL
******
`http://bio-bwa.sourceforge.net/bwa.shtml#3 <http://bio-bwa.sourceforge.net/bwa.shtml#3>`_

Docstring
*********
bwa - Burrows-Wheeler Alignment Tool

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

Outputs
-------
======  ======  ===============
name    type    documentation
======  ======  ===============
out     Stdout
======  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

=========  =============  ========  ==========  ===============
name       type           prefix      position  documentation
=========  =============  ========  ==========  ===============
reference  FastaWithDict                     9
reads      Fastq                            10
=========  =============  ========  ==========  ===============

Optional inputs
***************

===========================  ==================  ========  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                         type                prefix      position  documentation
===========================  ==================  ========  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
threads                      Optional<Integer>   -t                    Number of threads. (default = 1)
minimumSeedLength            Optional<Integer>   -k                    Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. (Default: 19)
bandwidth                    Optional<Integer>   -w                    Essentially, gaps longer than ${bandWidth} will not be found. Note that the maximum gap length is also affected by the scoring matrix and the hit length, not solely determined by this option. (Default: 100)
offDiagonalXDropoff          Optional<Integer>   -d                    (Z-dropoff): Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. (Default: 100)
reseedTrigger                Optional<Float>     -r                    Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. (Default: 1.5)
occurenceDiscard             Optional<Integer>   -c                    Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. (Default: 10000)
performSW                    Optional<Boolean>   -P                    In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.
matchingScore                Optional<Integer>   -A                    Matching score. (Default: 1)
mismatchPenalty              Optional<Integer>   -B                    Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. (Default: 4)
openGapPenalty               Optional<Integer>   -O                    Gap open penalty. (Default: 6)
gapExtensionPenalty          Optional<Integer>   -E                    Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). (Default: 1)
clippingPenalty              Optional<Integer>   -L                    Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. (Default: 5)
unpairedReadPenalty          Optional<Integer>   -U                    Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. It compares these two scores to determine whether we should force pairing. (Default: 9)
assumeInterleavedFirstInput  Optional<Boolean>   -p                    Assume the first input query file is interleaved paired-end FASTA/Q.
readGroupHeaderLine          Optional<String>    -R                    Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’. (Default=null)
outputAlignmentThreshold     Optional<Integer>   -T                    Don’t output alignment with score lower than INT. Only affects output. (Default: 30)
outputAllElements            Optional<Boolean>   -a                    Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.
appendComments               Optional<Boolean>   -C                    Append append FASTA/Q comment to SAM output. This option can be used to transfer read meta information (e.g. barcode) to the SAM output. Note that the FASTA/Q comment (the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). Malformated comments lead to incorrect SAM output.
hardClipping                 Optional<Boolean>   -H                    Use hard clipping ’H’ in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences.
markShorterSplits            Optional<Boolean>   -M                    Mark shorter split hits as secondary (for Picard compatibility).
verboseLevel                 Optional<Integer>   -v                    Control the verbose level of the output. This option has not been fully supported throughout BWA. Ideally, a value: 0 for disabling all the output to stderr; 1 for outputting errors only; 2 for warnings and errors; 3 for all normal messages; 4 or higher for debugging. When this option takes value 4, the output is not SAM. (Default: 3)
mates                        Optional<Fastq>                       11
outputFilename               Optional<Filename>
===========================  ==================  ========  ==========  =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: Michael Franklin


*BWA-MEM was last updated on 2018-12-24*.
*This page was automatically generated on 2019-05-03*.
