
.. include:: samtoolsview_v1.7.0
.. include:: samtoolsview_1.9

SamTools: View
=============================

Description
-------------

Tool identifier: ``SamToolsView``

Tool path: ``janis_bioinformatics.tools.samtools.view.latest import SamToolsViewLatest``

Version: v1.7.0

Docker: ``biocontainers/samtools:v1.7.0_cv3``

Versions
*********

- v1.7.0 (current)
- `1.9 <samtoolsview_1.9.html>`_

Documentation
-------------

URL
******
`http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS <http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS>`_

Tool documentation
******************
Ensure SAMTOOLS.SORT is inheriting from parent metadata
        
---------------------------------------------------------------------------------------------------
    
With no options or regions specified, prints all alignments in the specified input alignment file 
(in SAM, BAM, or CRAM format) to standard output in SAM format (with no header).

You may specify one or more space-separated region specifications after the input filename to 
restrict output to only those alignments which overlap the specified region(s). 
Use of region specifications requires a coordinate-sorted and indexed input file (in BAM or CRAM format).

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

======  ======  ========  ==========  ===============
name    type    prefix      position  documentation
======  ======  ========  ==========  ===============
sam     SAM                       10
======  ======  ========  ==========  ===============

Optional inputs
***************

=====================================  =======================  ========  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                                   type                     prefix      position  documentation
=====================================  =======================  ========  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================
cramOutput                             Optional<Boolean>        -C                 5  Output in the CRAM format (requires -T).
compressedBam                          Optional<Boolean>        -1                 5  Enable fast BAM compression (implies -b).
uncompressedBam                        Optional<Boolean>        -u                 5  Output uncompressed BAM. This option saves time spent on compression/decompression and is thus preferred when the output is piped to another samtools command.
onlyOutputHeader                       Optional<Boolean>        -H                 5  Output the header only.
countAlignments                        Optional<Boolean>        -c                 5  Instead of printing the alignments, only count them and print the total number. All filter options, such as -f, -F, and -q, are taken into account.
writeAlignments                        Optional<File>           -U                 5  Write alignments that are not selected by the various filter options to FILE. When this option is used, all alignments (or all alignments intersecting the regions specified) are written to either the output file or this file, but never both.
inputTSV                               Optional<File>           -t                 5  A tab-delimited FILE. Each line must contain the reference name in the first column and the length of the reference in the second column, with one line for each distinct reference. Any additional fields beyond the second column are ignored. This file also defines the order of the reference sequences in sorting. If you run: `samtools faidx <ref.fa>', the resulting index file <ref.fa>.fai can be used as this FILE.
onlyOverlapping                        Optional<File>           -L                 5  Only output alignments overlapping the input BED FILE [null].
useMultiRegionIterator                 Optional<Boolean>        -M                 5  Use the multi-region iterator on the union of the BED file and command-line region arguments. This avoids re-reading the same regions of files so can sometimes be much faster. Note this also removes duplicate sequences. Without this a sequence that overlaps multiple regions specified on the command line will be reported multiple times.
outputAlignmentsInReadGroup            Optional<String>         -r                 5  Output alignments in read group STR [null]. Note that records with no RG tag will also be output when using this option. This behaviour may change in a future release.
outputAlignmentsInFileReadGroups       Optional<File>           -R                 5  Output alignments in read groups listed in FILE [null]. Note that records with no RG tag will also be output when using this option. This behaviour may change in a future release.
mapqThreshold                          Optional<Integer>        -q                 5  Skip alignments with MAPQ smaller than INT [0].
outputAlignmentsInLibrary              Optional<String>         -l                 5  Only output alignments in library STR [null].
outputAlignmentsMeetingCIGARThreshold  Optional<Integer>        -m                 5  Only output alignments with number of CIGAR bases consuming query sequence ≥ INT [0]
outputAlignmentsWithBitsSet            Optional<String>         -f                 5  Only output alignments with all bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].
doNotOutputAlignmentsWithBitsSet       Optional<String>         -F                 5  Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].
doNotOutputAlignmentsWithAllBitsSet    Optional<String>         -G                 5  Do not output alignments with all bits set in INT present in the FLAG field. This is the opposite of -f such that -f12 -G12 is the same as no filtering at all. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].
readTagToExclude                       Optional<String>         -x                 5  Read tag to exclude from output (repeatable) [null]
collapseBackwardCIGAR                  Optional<Boolean>        -B                 5  Collapse the backward CIGAR operation.
subsamplingProportion                  Optional<Float>          -s                 5  Output only a proportion of the input alignments. This subsampling acts in the same way on all of the alignment records in the same template or read pair, so it never keeps a read but not its mate. The integer and fractional parts of the -s INT.FRAC option are used separately: the part after the decimal point sets the fraction of templates/pairs to be kept, while the integer part is used as a seed that influences which subset of reads is kept.
threads                                Optional<Integer>        -@                 5  Number of BAM compression threads to use in addition to main thread [0].
reference                              Optional<FastaWithDict>  -T                 6  A FASTA format reference FILE, optionally compressed by bgzip and ideally indexed by samtools faidx. If an index is not present, one will be generated for you.
outputFilename                         Optional<Filename>       -o                 5  Output to FILE [stdout].
=====================================  =======================  ========  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: Michael Franklin


*SamTools: View was last updated on 2019-01-24*.
*This page was automatically generated on 2019-07-24*.
