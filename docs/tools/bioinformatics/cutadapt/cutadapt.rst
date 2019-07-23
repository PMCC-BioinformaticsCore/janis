
.. include:: cutadapt_1.18

Cutadapt
===================
Tool identifier: ``cutadapt``
Tool path: ``janis_bioinformatics.tools.cutadapt.cutadapt_1_18 import CutAdapt_1_18``

Version: 1.18
Docker: ``quay.io/biocontainers/cutadapt:1.18--py37h14c3975_1``



Documentation
-------------

URL
******
`https://cutadapt.readthedocs.io/en/stable/ <https://cutadapt.readthedocs.io/en/stable/>`_

Description
*********

Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence 
from your high-throughput sequencing reads.

Cleaning your data in this way is often required: Reads from small-RNA sequencing 
contain the 3’ sequencing adapter because the read is longer than the molecule that is sequenced. 
Amplicon reads start with a primer sequence. Poly-A tails are useful for pulling out RNA from your sample, 
but often you don’t want them to be in your reads.
Cutadapt helps with these trimming tasks by finding the adapter or primer sequences in an error-tolerant way. 
It can also modify and filter reads in various ways. Adapter sequences can contain IUPAC wildcard characters. 
Also, paired-end reads and even colorspace data is supported. If you want, you can also just demultiplex your 
input data, without removing adapter sequences at all.

Cutadapt comes with an extensive suite of automated tests and is available under the terms of the MIT license.
If you use Cutadapt, please cite DOI:10.14806/ej.17.1.200 .

Outputs
-------
======  ======  ===============
name    type    documentation
======  ======  ===============
out     Fastq
======  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======  ======  ========  ==========  ===============
name    type    prefix      position  documentation
======  ======  ========  ==========  ===============
fastq   Fastq                      5
======  ======  ========  ==========  ===============

Optional inputs
***************

==========================  ==================  ==========================  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================
name                        type                prefix                      position    documentation
==========================  ==================  ==========================  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================
adapter                     Optional<String>    -a                                      Sequence of an adapter ligated to the 3' end (paired data: of the first read). The adapter and subsequent bases are trimmed. If a '$' character is appended ('anchoring'), the adapter is only found if it is a suffix of the read.
outputFilename              Optional<Filename>  -o                                      Write trimmed reads to FILE. FASTQ or FASTA format is chosen depending on input. The summary report is sent to standard output. Use '{name}' in FILE to demultiplex reads into multiple files. Default: write to standard output
secondReadFile              Optional<Filename>  -p                                      Write second read in a pair to FILE.
debug                       Optional<Boolean>   --debug                                 Print debugging information.
noIndels                    Optional<Boolean>   --no-indels                             Allow only mismatches in alignments. Default: allow both mismatches and indels
matchReadWildcards          Optional<Boolean>   --match-read-wildcards                  Interpret IUPAC wildcards in reads. Default: False
trimN                       Optional<Boolean>   --trim-n                                Trim N's on ends of reads.
discardCasava               Optional<Boolean>   --discard-casava                        Discard reads that did not pass CASAVA filtering (header has :Y:).
quiet                       Optional<Boolean>   --quiet                                 Print only error messages.
stripF3                     Optional<Boolean>   --strip-f3                              Strip the _F3 suffix of read names
noZeroCap                   Optional<Boolean>   --no-zero-cap                           Disable zero capping
interleaved                 Optional<Boolean>   --interleaved                           Read and write interleaved paired-end reads.
discardTrimmed              Optional<Boolean>   --discard-trimmed                       Discard reads that contain an adapter. Also use -O to avoid discarding too many randomly matching reads!
discardUntrimmed            Optional<Boolean>   --discard-untrimmed                     Discard reads that do not contain an adapter.
maq                         Optional<Boolean>   --maq                                   MAQ- and BWA-compatible colorspace output. This enables -c, -d, -t, --strip-f3 and -y '/1'.
pairFilter                  Optional<String>    --pair-filter=                          (any|both|first) Which of the reads in a paired-end read have to match the filtering criterion in order for the pair to be filtered. Default: any
nextseqTrim                 Optional<String>    --nextseq-trim=                         NextSeq-specific quality trimming (each read). Trims also dark cycles appearing as high-quality G bases.
action                      Optional<String>    --action=                               What to do with found adapters. trim: remove; mask: replace with 'N' characters; none: leave unchanged (useful with --discard-untrimmed). Default: trim
qualityBase                 Optional<String>    --quality-base=                         Assume that quality values in FASTQ are encoded as ascii(quality + N). This needs to be set to 64 for some old Illumina FASTQ files. Default: 33
lengthTag                   Optional<String>    --length-tag=                           Search for TAG followed by a decimal number in the description field of the read. Replace the decimal number with the correct length of the trimmed read. For example, use --length-tag 'length=' to correct fields like 'length=123'.
stripSuffix                 Optional<String>    --strip-suffix=                         Remove this suffix from read names if present. Can be given multiple times.
maxN                        Optional<Integer>   --max-n=                                Discard reads with more than COUNT 'N' bases. If COUNT is a number between 0 and 1, it is interpreted as a fraction of the read length.
report                      Optional<String>    --report=                               Which type of report to print. Default: full
infoFile                    Optional<String>    --info-file=                            Write information about each read and its adapter matches into FILE. See the documentation for the file format.
wildcardFile                Optional<String>    --wildcard-file=                        When the adapter has N wildcard bases, write adapter bases matching wildcard positions to FILE. (Inaccurate with indels.)
tooShortOutput              Optional<String>    --too-short-output=                     Write reads that are too short (according to length specified by -m) to FILE. Default: discard reads
tooLongOutput               Optional<String>    --too-long-output=                      Write reads that are too long (according to length specified by -M) to FILE. Default: discard reads
untrimmedOutput             Optional<String>    --untrimmed-output=                     Write reads that do not contain any adapter to FILE. Default: output to same file as trimmed reads
untrimmedPairedOutput       Optional<String>    --untrimmed-paired-output=              Write second read in a pair to this FILE when no adapter was found. Use with --untrimmed-output. Default: output to same file as trimmed reads
tooShortPairedOutput        Optional<String>    --too-short-paired-output=              Write second read in a pair to this file if pair is too short. Use also --too-short-output.
tooLongPairedOutput         Optional<String>    --too-long-paired-output=               Write second read in a pair to this file if pair is too long. Use also --too-long-output.
inputFileFormat             Optional<String>    -f                                      Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. Ignored when reading csfasta/qual files. Default: auto-detect from file name extension.
cores                       Optional<Integer>   -j                                      Number of CPU cores to use. Use 0 to auto-detect. Default: 1
adapter_g                   Optional<String>    -g                                      Sequence of an adapter ligated to the 5' end (paired data: of the first read). The adapter and any preceding bases are trimmed. Partial matches at the 5' end are allowed. If a '^' character is prepended ('anchoring'), the adapter is only found if it is a prefix of the read.
adapter_both                Optional<String>    -b                                      Sequence of an adapter that may be ligated to the 5' or 3' end (paired data: of the first read). Both types of matches as described under -a und -g are allowed. If the first base of the read is part of the match, the behavior is as with -g, otherwise as with -a. This option is mostly for rescuing failed library preparations - do not use if you know which end your adapter was ligated to!
maximumErrorRate            Optional<Float>     -e                                      Maximum allowed error rate as value between 0 and 1 (no. of errors divided by length of matching region). Default: 0.1 (=10%)
removeNAdapters             Optional<Integer>   -n                                      Remove up to COUNT adapters from each read. Default: 1
overlapRequirement          Optional<Integer>   -O                                      Require MINLENGTH overlap between read and adapter for an adapter to be found. Default: 3
removeNBases                Optional<Integer>   -u                                      Remove bases from each read (first read only if paired). If LENGTH is positive, remove bases from the beginning. If LENGTH is negative, remove bases from the end. Can be used twice if LENGTHs have different signs. This is applied *before* adapter trimming.
qualityCutoff               Optional<Integer>   -q                                      --quality-cutoff=[5'CUTOFF,]3'CUTOFF Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. Applied to both reads if data is paired. If one value is given, only the 3' end is trimmed. If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second.
shortenReadsToLength        Optional<Integer>   -l                                      Shorten reads to LENGTH. Positive values remove bases at the end while negative ones remove bases at the beginning. This and the following modifications are applied after adapter trimming.
readNamesPrefix             Optional<String>    -x                                      Add this prefix to read names. Use {name} to insert the name of the matching adapter.
readNamesSuffix             Optional<String>    -y                                      Add this suffix to read names; can also include {name}
minReadLength               Optional<Integer>   -m                                      --minimum-length=LEN[:LEN2] Discard reads shorter than LEN. Default: 0
maxReadsLength              Optional<Integer>   -M                                      --maximum-length=LEN[:LEN2] Discard reads longer than LEN. Default: no limit
middleReadMatchFile         Optional<String>    -r                                      When the adapter matches in the middle of a read, write the rest (after the adapter) to FILE.
removeMiddle3Adapter        Optional<String>    -A                                      3' adapter to be removed from second read in a pair.
removeMiddle5Adapter        Optional<String>    -G                                      5' adapter to be removed from second read in a pair.
removeMiddleBothAdapter     Optional<String>    -B                                      5'/3 adapter to be removed from second read in a pair.
removeNBasesFromSecondRead  Optional<Integer>   -U                                      Remove LENGTH bases from second read in a pair.
noMatchAdapterWildcards     Optional<Boolean>   -N                                      Do not interpret IUPAC wildcards in adapters.
colorspace                  Optional<Boolean>   -c                                      Enable colorspace mode
doubleEncode                Optional<Boolean>   -d                                      Double-encode colors (map 0,1,2,3,4 to A,C,G,T,N).
trimPrimer                  Optional<Boolean>   -t                                      Trim primer base and the first color
zeroCap                     Optional<Boolean>   -z                                      Change negative quality values to zero. Enabled by default in colorspace mode since many tools have problems with negative qualities
==========================  ==================  ==========================  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: **Unknown**


*Cutadapt was last updated on 2019-03-29*.
*This page was automatically generated on 2019-07-23*.
