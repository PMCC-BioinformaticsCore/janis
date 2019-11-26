:orphan:

Cutadapt
===================

1 contributor · 4 versions

:ID: ``cutadapt``
:Python: ``janis_bioinformatics.tools.cutadapt.versions import CutAdapt_1_18``
:Versions: 2.6, 2.5, 2.4, 1.18
:Container: quay.io/biocontainers/cutadapt:1.18--py37h14c3975_1
:Authors: Michael Franklin
:Citations: Martin, Marcel. “Cutadapt Removes Adapter Sequences from High-Throughput Sequencing Reads.” EMBnet.journal, vol. 17, no. 1, EMBnet Stichting, May 2011, p. 10. Crossref, doi:10.14806/ej.17.1.200.
:DOI: DOI:10.14806/ej.17.1.200
:Created: 2019-03-21
:Updated: 2019-03-29
:Required inputs:
   - ``fastq: FastqGzPair``
:Outputs: 
   - ``out: FastqGzPair``

Documentation
-------------

URL: `https://cutadapt.readthedocs.io/en/stable/ <https://cutadapt.readthedocs.io/en/stable/>`_


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

------

Additional configuration (inputs)
---------------------------------

==========================  ==================  =====================================================================================================================================================================================================================================================================================================================================================================================================
name                        type                documentation
==========================  ==================  =====================================================================================================================================================================================================================================================================================================================================================================================================
fastq                       FastqGzPair
adapter                     Optional<String>    Sequence of an adapter ligated to the 3' end (paired data: of the first read). The adapter and subsequent bases are trimmed. If a '$' character is appended ('anchoring'), the adapter is only found if it is a suffix of the read.
outputFilename              Optional<Filename>  Write trimmed reads to FILE. FASTQ or FASTA format is chosen depending on input. The summary report is sent to standard output. Use '{name}' in FILE to demultiplex reads into multiple files. Default: write to standard output
secondReadFile              Optional<Filename>  Write second read in a pair to FILE.
debug                       Optional<Boolean>   Print debugging information.
noIndels                    Optional<Boolean>   Allow only mismatches in alignments. Default: allow both mismatches and indels
matchReadWildcards          Optional<Boolean>   Interpret IUPAC wildcards in reads. Default: False
trimN                       Optional<Boolean>   Trim N's on ends of reads.
discardCasava               Optional<Boolean>   Discard reads that did not pass CASAVA filtering (header has :Y:).
quiet                       Optional<Boolean>   Print only error messages.
stripF3                     Optional<Boolean>   Strip the _F3 suffix of read names
noZeroCap                   Optional<Boolean>   Disable zero capping
interleaved                 Optional<Boolean>   Read and write interleaved paired-end reads.
discardTrimmed              Optional<Boolean>   Discard reads that contain an adapter. Also use -O to avoid discarding too many randomly matching reads!
discardUntrimmed            Optional<Boolean>   Discard reads that do not contain an adapter.
maq                         Optional<Boolean>   MAQ- and BWA-compatible colorspace output. This enables -c, -d, -t, --strip-f3 and -y '/1'.
pairFilter                  Optional<String>    (any|both|first) Which of the reads in a paired-end read have to match the filtering criterion in order for the pair to be filtered. Default: any
nextseqTrim                 Optional<String>    NextSeq-specific quality trimming (each read). Trims also dark cycles appearing as high-quality G bases.
action                      Optional<String>    What to do with found adapters. trim: remove; mask: replace with 'N' characters; none: leave unchanged (useful with --discard-untrimmed). Default: trim
qualityBase                 Optional<String>    Assume that quality values in FASTQ are encoded as ascii(quality + N). This needs to be set to 64 for some old Illumina FASTQ files. Default: 33
lengthTag                   Optional<String>    Search for TAG followed by a decimal number in the description field of the read. Replace the decimal number with the correct length of the trimmed read. For example, use --length-tag 'length=' to correct fields like 'length=123'.
stripSuffix                 Optional<String>    Remove this suffix from read names if present. Can be given multiple times.
maxN                        Optional<Integer>   Discard reads with more than COUNT 'N' bases. If COUNT is a number between 0 and 1, it is interpreted as a fraction of the read length.
report                      Optional<String>    Which type of report to print. Default: full
infoFile                    Optional<String>    Write information about each read and its adapter matches into FILE. See the documentation for the file format.
wildcardFile                Optional<String>    When the adapter has N wildcard bases, write adapter bases matching wildcard positions to FILE. (Inaccurate with indels.)
tooShortOutput              Optional<String>    Write reads that are too short (according to length specified by -m) to FILE. Default: discard reads
tooLongOutput               Optional<String>    Write reads that are too long (according to length specified by -M) to FILE. Default: discard reads
untrimmedOutput             Optional<String>    Write reads that do not contain any adapter to FILE. Default: output to same file as trimmed reads
untrimmedPairedOutput       Optional<String>    Write second read in a pair to this FILE when no adapter was found. Use with --untrimmed-output. Default: output to same file as trimmed reads
tooShortPairedOutput        Optional<String>    Write second read in a pair to this file if pair is too short. Use also --too-short-output.
tooLongPairedOutput         Optional<String>    Write second read in a pair to this file if pair is too long. Use also --too-long-output.
inputFileFormat             Optional<String>    Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. Ignored when reading csfasta/qual files. Default: auto-detect from file name extension.
cores                       Optional<Integer>   Number of CPU cores to use. Use 0 to auto-detect. Default: 1
adapter_g                   Optional<String>    Sequence of an adapter ligated to the 5' end (paired data: of the first read). The adapter and any preceding bases are trimmed. Partial matches at the 5' end are allowed. If a '^' character is prepended ('anchoring'), the adapter is only found if it is a prefix of the read.
adapter_both                Optional<String>    Sequence of an adapter that may be ligated to the 5' or 3' end (paired data: of the first read). Both types of matches as described under -a und -g are allowed. If the first base of the read is part of the match, the behavior is as with -g, otherwise as with -a. This option is mostly for rescuing failed library preparations - do not use if you know which end your adapter was ligated to!
maximumErrorRate            Optional<Float>     Maximum allowed error rate as value between 0 and 1 (no. of errors divided by length of matching region). Default: 0.1 (=10%)
removeNAdapters             Optional<Integer>   Remove up to COUNT adapters from each read. Default: 1
overlapRequirement          Optional<Integer>   Require MINLENGTH overlap between read and adapter for an adapter to be found. Default: 3
removeNBases                Optional<Integer>   Remove bases from each read (first read only if paired). If LENGTH is positive, remove bases from the beginning. If LENGTH is negative, remove bases from the end. Can be used twice if LENGTHs have different signs. This is applied *before* adapter trimming.
qualityCutoff               Optional<Integer>   --quality-cutoff=[5'CUTOFF,]3'CUTOFF Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. Applied to both reads if data is paired. If one value is given, only the 3' end is trimmed. If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second.
shortenReadsToLength        Optional<Integer>   Shorten reads to LENGTH. Positive values remove bases at the end while negative ones remove bases at the beginning. This and the following modifications are applied after adapter trimming.
readNamesPrefix             Optional<String>    Add this prefix to read names. Use {name} to insert the name of the matching adapter.
readNamesSuffix             Optional<String>    Add this suffix to read names; can also include {name}
minReadLength               Optional<Integer>   --minimum-length=LEN[:LEN2] Discard reads shorter than LEN. Default: 0
maxReadsLength              Optional<Integer>   --maximum-length=LEN[:LEN2] Discard reads longer than LEN. Default: no limit
middleReadMatchFile         Optional<String>    When the adapter matches in the middle of a read, write the rest (after the adapter) to FILE.
removeMiddle3Adapter        Optional<String>    3' adapter to be removed from second read in a pair.
removeMiddle5Adapter        Optional<String>    5' adapter to be removed from second read in a pair.
removeMiddleBothAdapter     Optional<String>    5'/3 adapter to be removed from second read in a pair.
removeNBasesFromSecondRead  Optional<Integer>   Remove LENGTH bases from second read in a pair.
noMatchAdapterWildcards     Optional<Boolean>   Do not interpret IUPAC wildcards in adapters.
colorspace                  Optional<Boolean>   Enable colorspace mode
doubleEncode                Optional<Boolean>   Double-encode colors (map 0,1,2,3,4 to A,C,G,T,N).
trimPrimer                  Optional<Boolean>   Trim primer base and the first color
zeroCap                     Optional<Boolean>   Change negative quality values to zero. Enabled by default in colorspace mode since many tools have problems with negative qualities
==========================  ==================  =====================================================================================================================================================================================================================================================================================================================================================================================================

