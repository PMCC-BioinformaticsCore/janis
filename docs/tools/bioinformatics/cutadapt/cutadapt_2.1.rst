:orphan:

Cutadapt
===================

``cutadapt`` · *1 contributor · 5 versions*

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

       from janis_bioinformatics.tools.cutadapt.versions import CutAdapt_2_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "cutadapt_step",
           CutAdapt_2_1(
               outputPrefix=None,
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
       outputPrefix: <value>




5. Run cutadapt with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       cutadapt





Information
------------

:ID: ``cutadapt``
:URL: `https://cutadapt.readthedocs.io/en/stable/ <https://cutadapt.readthedocs.io/en/stable/>`_
:Versions: 2.6, 2.5, 2.4, 2.1, 1.18
:Container: quay.io/biocontainers/cutadapt:2.1--py37h14c3975_0
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
outputPrefix                String                                                             Used for naming purposes
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

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task cutadapt {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       String outputPrefix
       Array[File] fastq
       Array[String]? adapter
       String? outputFilename
       String? secondReadFile
       Int? cores
       String? front
       String? anywhere
       Float? errorRate
       Boolean? noIndels
       Int? times
       Int? overlap
       Boolean? matchReadWildcards
       Boolean? noMatchAdapterWildcards
       String? action
       Int? cut
       String? nextseqTrim
       Int? qualityCutoff
       Boolean? qualityBase
       Int? length
       Int? trimN
       Int? lengthTag
       String? stripSuffix
       String? prefix
       String? suffix
       Boolean? zeroCap
       Int? minimumLength
       Int? maximumLength
       Float? maxN
       Boolean? discardTrimmed
       Boolean? discardUntrimmed
       Boolean? discardCasava
       Boolean? quiet
       String? compressionLevel
       String? infoFile
       String? restFile
       String? wildcardFile
       String? tooShortOutput
       String? tooLongOutput
       String? untrimmedOutput
       Array[String]? removeMiddle3Adapter
       String? removeMiddle5Adapter
       String? removeMiddleBothAdapter
       String? removeNBasesFromSecondRead
       String? pairAdapters
       String? pairFilter
       Boolean? interleaved
       String? untrimmedPairedOutput
       String? tooShortPairedOutput
       String? tooLongPairedOutput
     }
     command <<<
       set -e
       cutadapt \
         ~{if (defined(adapter) && length(select_first([adapter])) > 0) then "-a '" + sep("' -a '", select_first([adapter])) + "'" else ""} \
         -o '~{select_first([outputFilename, "~{outputPrefix}-R1.fastq.gz"])}' \
         -p '~{select_first([secondReadFile, "~{outputPrefix}-R2.fastq.gz"])}' \
         ~{if defined(cores) then ("--cores " + cores) else ''} \
         ~{if defined(front) then ("--front '" + front + "'") else ""} \
         ~{if defined(anywhere) then ("--anywhere '" + anywhere + "'") else ""} \
         ~{if defined(errorRate) then ("--error-rate " + errorRate) else ''} \
         ~{if (defined(noIndels) && select_first([noIndels])) then "--no-indels" else ""} \
         ~{if defined(times) then ("--times " + times) else ''} \
         ~{if defined(overlap) then ("--overlap " + overlap) else ''} \
         ~{if (defined(matchReadWildcards) && select_first([matchReadWildcards])) then "--match-read-wildcards" else ""} \
         ~{if (defined(noMatchAdapterWildcards) && select_first([noMatchAdapterWildcards])) then "--no-match-adapter-wildcards" else ""} \
         ~{if defined(action) then ("--action '" + action + "'") else ""} \
         ~{if defined(cut) then ("--cut " + cut) else ''} \
         ~{if defined(nextseqTrim) then ("--nextseq-trim '" + nextseqTrim + "'") else ""} \
         ~{if defined(qualityCutoff) then ("--quality-cutoff " + qualityCutoff) else ''} \
         ~{if (defined(qualityBase) && select_first([qualityBase])) then "--quality-base" else ""} \
         ~{if defined(length) then ("--length " + length) else ''} \
         ~{if defined(trimN) then ("--trim-n " + trimN) else ''} \
         ~{if defined(lengthTag) then ("--length-tag " + lengthTag) else ''} \
         ~{if defined(stripSuffix) then ("--strip-suffix '" + stripSuffix + "'") else ""} \
         ~{if defined(prefix) then ("--prefix '" + prefix + "'") else ""} \
         ~{if defined(suffix) then ("--suffix '" + suffix + "'") else ""} \
         ~{if (defined(zeroCap) && select_first([zeroCap])) then "--zero-cap" else ""} \
         ~{if defined(minimumLength) then ("--minimum-length " + minimumLength) else ''} \
         ~{if defined(maximumLength) then ("--maximum-length " + maximumLength) else ''} \
         ~{if defined(maxN) then ("--max-n " + maxN) else ''} \
         ~{if (defined(discardTrimmed) && select_first([discardTrimmed])) then "--discard-trimmed" else ""} \
         ~{if (defined(discardUntrimmed) && select_first([discardUntrimmed])) then "--discard-untrimmed" else ""} \
         ~{if (defined(discardCasava) && select_first([discardCasava])) then "--discard-casava" else ""} \
         ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
         ~{if defined(compressionLevel) then ("-Z '" + compressionLevel + "'") else ""} \
         ~{if defined(infoFile) then ("--info-file '" + infoFile + "'") else ""} \
         ~{if defined(restFile) then ("--rest-file '" + restFile + "'") else ""} \
         ~{if defined(wildcardFile) then ("--wildcard-file '" + wildcardFile + "'") else ""} \
         ~{if defined(tooShortOutput) then ("--too-short-output '" + tooShortOutput + "'") else ""} \
         ~{if defined(tooLongOutput) then ("--too-long-output '" + tooLongOutput + "'") else ""} \
         ~{if defined(untrimmedOutput) then ("--untrimmed-output '" + untrimmedOutput + "'") else ""} \
         ~{if (defined(removeMiddle3Adapter) && length(select_first([removeMiddle3Adapter])) > 0) then "-A '" + sep("' -A '", select_first([removeMiddle3Adapter])) + "'" else ""} \
         ~{if defined(removeMiddle5Adapter) then ("-G '" + removeMiddle5Adapter + "'") else ""} \
         ~{if defined(removeMiddleBothAdapter) then ("-B '" + removeMiddleBothAdapter + "'") else ""} \
         ~{if defined(removeNBasesFromSecondRead) then ("-U '" + removeNBasesFromSecondRead + "'") else ""} \
         ~{if defined(pairAdapters) then ("--pair-adapters '" + pairAdapters + "'") else ""} \
         ~{if defined(pairFilter) then ("--pair-filter '" + pairFilter + "'") else ""} \
         ~{if (defined(interleaved) && select_first([interleaved])) then "--interleaved" else ""} \
         ~{if defined(untrimmedPairedOutput) then ("--untrimmed-paired-output '" + untrimmedPairedOutput + "'") else ""} \
         ~{if defined(tooShortPairedOutput) then ("--too-short-paired-output '" + tooShortPairedOutput + "'") else ""} \
         ~{if defined(tooLongPairedOutput) then ("--too-long-paired-output '" + tooLongPairedOutput + "'") else ""} \
         ~{if length(fastq) > 0 then "'" + sep("' '", fastq) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 5, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/cutadapt:2.1--py37h14c3975_0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4, 4])}G"
       preemptible: 2
     }
     output {
       Array[File] out = [(outputPrefix + "-R1.fastq.gz"), (outputPrefix + "-R2.fastq.gz")]
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Cutadapt
   doc: |
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

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/cutadapt:2.1--py37h14c3975_0

   inputs:
   - id: outputPrefix
     label: outputPrefix
     doc: Used for naming purposes
     type: string
   - id: fastq
     label: fastq
     type:
       type: array
       items: File
     inputBinding:
       position: 5
   - id: adapter
     label: adapter
     doc: |-
       Sequence of an adapter ligated to the 3' end (paired data: of the first read). The adapter and subsequent bases are trimmed. If a '$' character is appended ('anchoring'), the adapter is only found if it is a suffix of the read.
     type:
     - type: array
       inputBinding:
         prefix: -a
       items: string
     - 'null'
     inputBinding: {}
   - id: outputFilename
     label: outputFilename
     doc: |-
       Write trimmed reads to FILE. FASTQ or FASTA format is chosen depending on input. The summary report is sent to standard output. Use '{name}' in FILE to demultiplex reads into multiple files. Default: write to standard output
     type:
     - string
     - 'null'
     default: generated-R1.fastq.gz
     inputBinding:
       prefix: -o
       valueFrom: $(inputs.outputPrefix)-R1.fastq.gz
   - id: secondReadFile
     label: secondReadFile
     doc: Write second read in a pair to FILE.
     type:
     - string
     - 'null'
     default: generated-R2.fastq.gz
     inputBinding:
       prefix: -p
       valueFrom: $(inputs.outputPrefix)-R2.fastq.gz
   - id: cores
     label: cores
     doc: '(-j)  Number of CPU cores to use. Use 0 to auto-detect. Default: 1'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --cores
       separate: true
   - id: front
     label: front
     doc: |-
       (-g)  Sequence of an adapter ligated to the 5' end (paired data: of the first read). The adapter and any preceding bases are trimmed. Partial matches at the 5' end are allowed. If a '^' character is prepended ('anchoring'), the adapter is only found if it is a prefix of the read.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --front
       separate: true
   - id: anywhere
     label: anywhere
     doc: |-
       (-b)  Sequence of an adapter that may be ligated to the 5' or 3' end (paired data: of the first read). Both types of matches as described under -a und -g are allowed. If the first base of the read is part of the match, the behavior is as with -g, otherwise as with -a. This option is mostly for rescuing failed library preparations - do not use if you know which end your adapter was ligated to!
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --anywhere
       separate: true
   - id: errorRate
     label: errorRate
     doc: |-
       (-e)  Maximum allowed error rate as value between 0 and 1 (no. of errors divided by length of matching region). Default: 0.1 (=10%)
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --error-rate
       separate: true
   - id: noIndels
     label: noIndels
     doc: 'Allow only mismatches in alignments. Default: allow both mismatches and indels'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-indels
       separate: true
   - id: times
     label: times
     doc: '(-n)  Remove up to COUNT adapters from each read. Default: 1'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --times
       separate: true
   - id: overlap
     label: overlap
     doc: |-
       (-O)  Require MINLENGTH overlap between read and adapter for an adapter to be found. Default: 3
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --overlap
       separate: true
   - id: matchReadWildcards
     label: matchReadWildcards
     doc: ' Interpret IUPAC wildcards in reads. Default: False'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --match-read-wildcards
       separate: true
   - id: noMatchAdapterWildcards
     label: noMatchAdapterWildcards
     doc: (-N)  Do not interpret IUPAC wildcards in adapters.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-match-adapter-wildcards
       separate: true
   - id: action
     label: action
     doc: |-
       (trim,mask,lowercase,none}  What to do with found adapters. mask: replace with 'N' characters; lowercase: convert to lowercase; none: leave unchanged (useful with --discard-untrimmed). Default: trim
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --action
       separate: true
   - id: cut
     label: cut
     doc: |-
       (-u)  Remove bases from each read (first read only if paired). If LENGTH is positive, remove bases from the beginning. If LENGTH is negative, remove bases from the end. Can be used twice if LENGTHs have different signs. This is applied *before* adapter trimming.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --cut
       separate: true
   - id: nextseqTrim
     label: nextseqTrim
     doc: |2-
        NextSeq-specific quality trimming (each read). Trims also dark cycles appearing as high-quality G bases.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --nextseq-trim
       separate: true
   - id: qualityCutoff
     label: qualityCutoff
     doc: |-
       (]3'CUTOFF, ]3'CUTOFF, -q)  Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. Applied to both reads if data is paired. If one value is given, only the 3' end is trimmed. If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --quality-cutoff
       separate: true
   - id: qualityBase
     label: qualityBase
     doc: |-
       Assume that quality values in FASTQ are encoded as ascii(quality + N). This needs to be set to 64 for some old Illumina FASTQ files. Default: 33
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --quality-base
       separate: true
   - id: length
     label: length
     doc: |-
       (-l)  Shorten reads to LENGTH. Positive values remove bases at the end while negative ones remove bases at the beginning. This and the following modifications are applied after adapter trimming.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --length
       separate: true
   - id: trimN
     label: trimN
     doc: Trim N's on ends of reads.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --trim-n
       separate: true
   - id: lengthTag
     label: lengthTag
     doc: |-
       Search for TAG followed by a decimal number in the description field of the read. Replace the decimal number with the correct length of the trimmed read. For example, use --length-tag 'length=' to correct fields like 'length=123'.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --length-tag
       separate: true
   - id: stripSuffix
     label: stripSuffix
     doc: ' Remove this suffix from read names if present. Can be given multiple times.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --strip-suffix
       separate: true
   - id: prefix
     label: prefix
     doc: |-
       (-x)  Add this prefix to read names. Use {name} to insert the name of the matching adapter.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --prefix
       separate: true
   - id: suffix
     label: suffix
     doc: (-y)  Add this suffix to read names; can also include {name}
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --suffix
       separate: true
   - id: zeroCap
     label: zeroCap
     doc: (-z) Change negative quality values to zero.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --zero-cap
       separate: true
   - id: minimumLength
     label: minimumLength
     doc: '(-m)  Discard reads shorter than LEN. Default: 0'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --minimum-length
       separate: true
   - id: maximumLength
     label: maximumLength
     doc: '(-M)  Discard reads longer than LEN. Default: no limit'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --maximum-length
       separate: true
   - id: maxN
     label: maxN
     doc: |-
       Discard reads with more than COUNT 'N' bases. If COUNT is a number between 0 and 1, it is interpreted as a fraction of the read length.
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --max-n
       separate: true
   - id: discardTrimmed
     label: discardTrimmed
     doc: |-
       (--discard)  Discard reads that contain an adapter. Use also -O to avoid discarding too many randomly matching reads.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --discard-trimmed
       separate: true
   - id: discardUntrimmed
     label: discardUntrimmed
     doc: (--trimmed-only)  Discard reads that do not contain an adapter.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --discard-untrimmed
       separate: true
   - id: discardCasava
     label: discardCasava
     doc: Discard reads that did not pass CASAVA filtering (header has :Y:).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --discard-casava
       separate: true
   - id: quiet
     label: quiet
     doc: |-
       Print only error messages. Which type of report to print: 'full' or 'minimal'. Default: full
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --quiet
       separate: true
   - id: compressionLevel
     label: compressionLevel
     doc: Use compression level 1 for gzipped output files (faster, but uses more space)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -Z
       separate: true
   - id: infoFile
     label: infoFile
     doc: |-
       Write information about each read and its adapter matches into FILE. See the documentation for the file format.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --info-file
       separate: true
   - id: restFile
     label: restFile
     doc: |-
       (-r)  When the adapter matches in the middle of a read, write the rest (after the adapter) to FILE.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --rest-file
       separate: true
   - id: wildcardFile
     label: wildcardFile
     doc: |-
       When the adapter has N wildcard bases, write adapter bases matching wildcard positions to FILE. (Inaccurate with indels.)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --wildcard-file
       separate: true
   - id: tooShortOutput
     label: tooShortOutput
     doc: |2-
        Write reads that are too short (according to length specified by -m) to FILE. Default: discard reads
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --too-short-output
       separate: true
   - id: tooLongOutput
     label: tooLongOutput
     doc: |2-
        Write reads that are too long (according to length specified by -M) to FILE. Default: discard reads
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --too-long-output
       separate: true
   - id: untrimmedOutput
     label: untrimmedOutput
     doc: |2-
        Write reads that do not contain any adapter to FILE. Default: output to same file as trimmed reads
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --untrimmed-output
       separate: true
   - id: removeMiddle3Adapter
     label: removeMiddle3Adapter
     doc: 3' adapter to be removed from second read in a pair.
     type:
     - type: array
       inputBinding:
         prefix: -A
         separate: true
       items: string
     - 'null'
     inputBinding: {}
   - id: removeMiddle5Adapter
     label: removeMiddle5Adapter
     doc: 5' adapter to be removed from second read in a pair.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -G
       separate: true
   - id: removeMiddleBothAdapter
     label: removeMiddleBothAdapter
     doc: 5'/3 adapter to be removed from second read in a pair.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -B
       separate: true
   - id: removeNBasesFromSecondRead
     label: removeNBasesFromSecondRead
     doc: Remove LENGTH bases from second read in a pair.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -U
       separate: true
   - id: pairAdapters
     label: pairAdapters
     doc: |-
       Treat adapters given with -a/-A etc. as pairs. Either both or none are removed from each read pair.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --pair-adapters
       separate: true
   - id: pairFilter
     label: pairFilter
     doc: |-
       {any,both,first} Which of the reads in a paired-end read have to match the filtering criterion in order for the pair to be filtered. Default: any
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --pair-filter
       separate: true
   - id: interleaved
     label: interleaved
     doc: Read and write interleaved paired-end reads.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --interleaved
       separate: true
   - id: untrimmedPairedOutput
     label: untrimmedPairedOutput
     doc: |2-
        Write second read in a pair to this FILE when no adapter was found. Use with --untrimmed-output. Default: output to same file as trimmed reads
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --untrimmed-paired-output
       separate: true
   - id: tooShortPairedOutput
     label: tooShortPairedOutput
     doc: |2-
        Write second read in a pair to this file if pair is too short. Use also --too-short-output.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --too-short-paired-output
       separate: true
   - id: tooLongPairedOutput
     label: tooLongPairedOutput
     doc: |2-
        Write second read in a pair to this file if pair is too long. Use also --too-long-output.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --too-long-paired-output
       separate: true

   outputs:
   - id: out
     label: out
     type:
       type: array
       items: File
     outputBinding:
       glob:
       - $((inputs.outputPrefix + "-R1.fastq.gz"))
       - $((inputs.outputPrefix + "-R2.fastq.gz"))
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: cutadapt
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: cutadapt


