:orphan:

Cutadapt
===================

*1 contributor · 4 versions*

cutadapt version 2.4
Copyright (C) 2010-2019 Marcel Martin <marcel.martin@scilifelab.se>
cutadapt removes adapter sequences from high-throughput sequencing reads.
Usage:
    cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq
For paired-end reads:
    cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq
Replace "ADAPTER" with the actual sequence of your 3' adapter. IUPAC wildcard
characters are supported. The reverse complement is *not* automatically
searched. All reads from input.fastq will be written to output.fastq with the
adapter sequence removed. Adapter matching is error-tolerant. Multiple adapter
sequences can be given (use further -a options), but only the best-matching
adapter will be removed.
Input may also be in FASTA format. Compressed input and output is supported and
auto-detected from the file name (.gz, .xz, .bz2). Use the file name '-' for
standard input/output. Without the -o option, output is sent to standard output.
Citation:
Marcel Martin. Cutadapt removes adapter sequences from high-throughput
sequencing reads. EMBnet.Journal, 17(1):10-12, May 2011.
http://dx.doi.org/10.14806/ej.17.1.200
Run "cutadapt - -help" to see all command-line options.
See https://cutadapt.readthedocs.io/ for full documentation.



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.cutadapt.versions import CutAdapt_2_5

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "cutadapt_step",
           CutAdapt_2_5(
               fastq=None,
           )
       )
       wf.output("out", source=cutadapt_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for cutadapt:

.. code-block:: bash

   # user inputs
   janis inputs cutadapt > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fastq:
       - fastq_0.fastq.gz
       - fastq_1.fastq.gz




5. Run cutadapt with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       cutadapt





Information
------------


:ID: ``cutadapt``
:URL: `https://cutadapt.readthedocs.io/en/stable/ <https://cutadapt.readthedocs.io/en/stable/>`_
:Versions: 2.6, 2.5, 2.4, 1.18
:Container: quay.io/biocontainers/cutadapt:2.5--py37h516909a_0
:Authors: Michael Franklin
:Citations: Martin, Marcel. “Cutadapt Removes Adapter Sequences from High-Throughput Sequencing Reads.” EMBnet.journal, vol. 17, no. 1, EMBnet Stichting, May 2011, p. 10. Crossref, doi:10.14806/ej.17.1.200.
:DOI: DOI:10.14806/ej.17.1.200
:Created: 2019-03-21
:Updated: 2019-07-23



Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     FastqGzPair
======  ===========  ===============



Additional configuration (inputs)
---------------------------------

==========================  =======================  ============================  ==========  ===========================================================================================================================================================================================================================================================================================================================================================================================================
name                        type                     prefix                          position  documentation
==========================  =======================  ============================  ==========  ===========================================================================================================================================================================================================================================================================================================================================================================================================
fastq                       FastqGzPair                                                     5
adapter                     Optional<Array<String>>  -a                                        Sequence of an adapter ligated to the 3' end (paired data: of the first read). The adapter and subsequent bases are trimmed. If a '$' character is appended ('anchoring'), the adapter is only found if it is a suffix of the read.
outputFilename              Optional<Filename>       -o                                        Write trimmed reads to FILE. FASTQ or FASTA format is chosen depending on input. The summary report is sent to standard output. Use '{name}' in FILE to demultiplex reads into multiple files. Default: write to standard output
secondReadFile              Optional<Filename>       -p                                        Write second read in a pair to FILE.
cores                       Optional<Integer>        --cores                                   (-j)  Number of CPU cores to use. Use 0 to auto-detect. Default: 1
front                       Optional<String>         --front                                   (-g)  Sequence of an adapter ligated to the 5' end (paired data: of the first read). The adapter and any preceding bases are trimmed. Partial matches at the 5' end are allowed. If a '^' character is prepended ('anchoring'), the adapter is only found if it is a prefix of the read.
anywhere                    Optional<String>         --anywhere                                (-b)  Sequence of an adapter that may be ligated to the 5' or 3' end (paired data: of the first read). Both types of matches as described under -a und -g are allowed. If the first base of the read is part of the match, the behavior is as with -g, otherwise as with -a. This option is mostly for rescuing failed library preparations - do not use if you know which end your adapter was ligated to!
errorRate                   Optional<Float>          --error-rate                              (-e)  Maximum allowed error rate as value between 0 and 1 (no. of errors divided by length of matching region). Default: 0.1 (=10%)
noIndels                    Optional<Boolean>        --no-indels                               Allow only mismatches in alignments. Default: allow both mismatches and indels
times                       Optional<Integer>        --times                                   (-n)  Remove up to COUNT adapters from each read. Default: 1
overlap                     Optional<Integer>        --overlap                                 (-O)  Require MINLENGTH overlap between read and adapter for an adapter to be found. Default: 3
matchReadWildcards          Optional<Boolean>        --match-read-wildcards                    Interpret IUPAC wildcards in reads. Default: False
noMatchAdapterWildcards     Optional<Boolean>        --no-match-adapter-wildcards              (-N)  Do not interpret IUPAC wildcards in adapters.
action                      Optional<String>         --action                                  (trim,mask,lowercase,none}  What to do with found adapters. mask: replace with 'N' characters; lowercase: convert to lowercase; none: leave unchanged (useful with --discard-untrimmed). Default: trim
cut                         Optional<Integer>        --cut                                     (-u)  Remove bases from each read (first read only if paired). If LENGTH is positive, remove bases from the beginning. If LENGTH is negative, remove bases from the end. Can be used twice if LENGTHs have different signs. This is applied *before* adapter trimming.
nextseqTrim                 Optional<String>         --nextseq-trim                            NextSeq-specific quality trimming (each read). Trims also dark cycles appearing as high-quality G bases.
qualityCutoff               Optional<Integer>        --quality-cutoff                          (]3'CUTOFF, ]3'CUTOFF, -q)  Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. Applied to both reads if data is paired. If one value is given, only the 3' end is trimmed. If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second.
qualityBase                 Optional<Boolean>        --quality-base                            Assume that quality values in FASTQ are encoded as ascii(quality + N). This needs to be set to 64 for some old Illumina FASTQ files. Default: 33
length                      Optional<Integer>        --length                                  (-l)  Shorten reads to LENGTH. Positive values remove bases at the end while negative ones remove bases at the beginning. This and the following modifications are applied after adapter trimming.
trimN                       Optional<Integer>        --trim-n                                  Trim N's on ends of reads.
lengthTag                   Optional<Integer>        --length-tag                              Search for TAG followed by a decimal number in the description field of the read. Replace the decimal number with the correct length of the trimmed read. For example, use --length-tag 'length=' to correct fields like 'length=123'.
stripSuffix                 Optional<String>         --strip-suffix                            Remove this suffix from read names if present. Can be given multiple times.
prefix                      Optional<String>         --prefix                                  (-x)  Add this prefix to read names. Use {name} to insert the name of the matching adapter.
suffix                      Optional<String>         --suffix                                  (-y)  Add this suffix to read names; can also include {name}
zeroCap                     Optional<Boolean>        --zero-cap                                (-z) Change negative quality values to zero.
minimumLength               Optional<Integer>        --minimum-length                          (-m)  Discard reads shorter than LEN. Default: 0
maximumLength               Optional<Integer>        --maximum-length                          (-M)  Discard reads longer than LEN. Default: no limit
maxN                        Optional<Float>          --max-n                                   Discard reads with more than COUNT 'N' bases. If COUNT is a number between 0 and 1, it is interpreted as a fraction of the read length.
discardTrimmed              Optional<Boolean>        --discard-trimmed                         (--discard)  Discard reads that contain an adapter. Use also -O to avoid discarding too many randomly matching reads.
discardUntrimmed            Optional<Boolean>        --discard-untrimmed                       (--trimmed-only)  Discard reads that do not contain an adapter.
discardCasava               Optional<Boolean>        --discard-casava                          Discard reads that did not pass CASAVA filtering (header has :Y:).
quiet                       Optional<Boolean>        --quiet                                   Print only error messages. Which type of report to print: 'full' or 'minimal'. Default: full
compressionLevel            Optional<String>         -Z                                        Use compression level 1 for gzipped output files (faster, but uses more space)
infoFile                    Optional<String>         --info-file                               Write information about each read and its adapter matches into FILE. See the documentation for the file format.
restFile                    Optional<String>         --rest-file                               (-r)  When the adapter matches in the middle of a read, write the rest (after the adapter) to FILE.
wildcardFile                Optional<String>         --wildcard-file                           When the adapter has N wildcard bases, write adapter bases matching wildcard positions to FILE. (Inaccurate with indels.)
tooShortOutput              Optional<String>         --too-short-output                        Write reads that are too short (according to length specified by -m) to FILE. Default: discard reads
tooLongOutput               Optional<String>         --too-long-output                         Write reads that are too long (according to length specified by -M) to FILE. Default: discard reads
untrimmedOutput             Optional<String>         --untrimmed-output                        Write reads that do not contain any adapter to FILE. Default: output to same file as trimmed reads
removeMiddle3Adapter        Optional<Array<String>>  -A                                        3' adapter to be removed from second read in a pair.
removeMiddle5Adapter        Optional<String>         -G                                        5' adapter to be removed from second read in a pair.
removeMiddleBothAdapter     Optional<String>         -B                                        5'/3 adapter to be removed from second read in a pair.
removeNBasesFromSecondRead  Optional<String>         -U                                        Remove LENGTH bases from second read in a pair.
pairAdapters                Optional<String>         --pair-adapters                           Treat adapters given with -a/-A etc. as pairs. Either both or none are removed from each read pair.
pairFilter                  Optional<String>         --pair-filter                             {any,both,first} Which of the reads in a paired-end read have to match the filtering criterion in order for the pair to be filtered. Default: any
interleaved                 Optional<Boolean>        --interleaved                             Read and write interleaved paired-end reads.
untrimmedPairedOutput       Optional<String>         --untrimmed-paired-output                 Write second read in a pair to this FILE when no adapter was found. Use with --untrimmed-output. Default: output to same file as trimmed reads
tooShortPairedOutput        Optional<String>         --too-short-paired-output                 Write second read in a pair to this file if pair is too short. Use also --too-short-output.
tooLongPairedOutput         Optional<String>         --too-long-paired-output                  Write second read in a pair to this file if pair is too long. Use also --too-long-output.
==========================  =======================  ============================  ==========  ===========================================================================================================================================================================================================================================================================================================================================================================================================
