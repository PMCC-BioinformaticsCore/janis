
bwamem
======
*bioinformatics*

Documentation
-------------

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

*This page was automatically generated*
