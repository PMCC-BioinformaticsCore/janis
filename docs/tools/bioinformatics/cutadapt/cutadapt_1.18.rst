:orphan:

Cutadapt
===================

``cutadapt`` · *1 contributor · 5 versions*


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


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.cutadapt.versions import CutAdapt_1_18

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "cutadapt_step",
           CutAdapt_1_18(
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
:Versions: 2.6, 2.5, 2.4, 2.1, 1.18
:Container: quay.io/biocontainers/cutadapt:1.18--py37h14c3975_1
:Authors: Michael Franklin
:Citations: Martin, Marcel. “Cutadapt Removes Adapter Sequences from High-Throughput Sequencing Reads.” EMBnet.journal, vol. 17, no. 1, EMBnet Stichting, May 2011, p. 10. Crossref, doi:10.14806/ej.17.1.200.
:DOI: DOI:10.14806/ej.17.1.200
:Created: 2019-03-21
:Updated: 2019-03-29


Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     FastqGzPair
======  ===========  ===============


Additional configuration (inputs)
---------------------------------

==========================  ==================  ==========================  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================
name                        type                prefix                        position  documentation
==========================  ==================  ==========================  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================
fastq                       FastqGzPair                                              5
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
       Array[File] fastq
       String? adapter
       String? outputFilename
       String? secondReadFile
       Boolean? debug
       Boolean? noIndels
       Boolean? matchReadWildcards
       Boolean? trimN
       Boolean? discardCasava
       Boolean? quiet
       Boolean? stripF3
       Boolean? noZeroCap
       Boolean? interleaved
       Boolean? discardTrimmed
       Boolean? discardUntrimmed
       Boolean? maq
       String? pairFilter
       String? nextseqTrim
       String? action
       String? qualityBase
       String? lengthTag
       String? stripSuffix
       Int? maxN
       String? report
       String? infoFile
       String? wildcardFile
       String? tooShortOutput
       String? tooLongOutput
       String? untrimmedOutput
       String? untrimmedPairedOutput
       String? tooShortPairedOutput
       String? tooLongPairedOutput
       String? inputFileFormat
       Int? cores
       String? adapter_g
       String? adapter_both
       Float? maximumErrorRate
       Int? removeNAdapters
       Int? overlapRequirement
       Int? removeNBases
       Int? qualityCutoff
       Int? shortenReadsToLength
       String? readNamesPrefix
       String? readNamesSuffix
       Int? minReadLength
       Int? maxReadsLength
       String? middleReadMatchFile
       String? removeMiddle3Adapter
       String? removeMiddle5Adapter
       String? removeMiddleBothAdapter
       Int? removeNBasesFromSecondRead
       Boolean? noMatchAdapterWildcards
       Boolean? colorspace
       Boolean? doubleEncode
       Boolean? trimPrimer
       Boolean? zeroCap
     }
     command <<<
       set -e
       cutadapt \
         ~{if defined(adapter) then ("-a '" + adapter + "'") else ""} \
         -o '~{select_first([outputFilename, "generated--R1.fastq.gz"])}' \
         -p '~{select_first([secondReadFile, "generated--R2.fastq.gz"])}' \
         ~{if (defined(debug) && select_first([debug])) then "--debug" else ""} \
         ~{if (defined(noIndels) && select_first([noIndels])) then "--no-indels" else ""} \
         ~{if (defined(matchReadWildcards) && select_first([matchReadWildcards])) then "--match-read-wildcards" else ""} \
         ~{if (defined(trimN) && select_first([trimN])) then "--trim-n" else ""} \
         ~{if (defined(discardCasava) && select_first([discardCasava])) then "--discard-casava" else ""} \
         ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
         ~{if (defined(stripF3) && select_first([stripF3])) then "--strip-f3" else ""} \
         ~{if (defined(noZeroCap) && select_first([noZeroCap])) then "--no-zero-cap" else ""} \
         ~{if (defined(interleaved) && select_first([interleaved])) then "--interleaved" else ""} \
         ~{if (defined(discardTrimmed) && select_first([discardTrimmed])) then "--discard-trimmed" else ""} \
         ~{if (defined(discardUntrimmed) && select_first([discardUntrimmed])) then "--discard-untrimmed" else ""} \
         ~{if (defined(maq) && select_first([maq])) then "--maq" else ""} \
         ~{if defined(pairFilter) then ("--pair-filter= '" + pairFilter + "'") else ""} \
         ~{if defined(nextseqTrim) then ("--nextseq-trim= '" + nextseqTrim + "'") else ""} \
         ~{if defined(action) then ("--action= '" + action + "'") else ""} \
         ~{if defined(qualityBase) then ("--quality-base= '" + qualityBase + "'") else ""} \
         ~{if defined(lengthTag) then ("--length-tag= '" + lengthTag + "'") else ""} \
         ~{if defined(stripSuffix) then ("--strip-suffix= '" + stripSuffix + "'") else ""} \
         ~{if defined(maxN) then ("--max-n= " + maxN) else ''} \
         ~{if defined(report) then ("--report= '" + report + "'") else ""} \
         ~{if defined(infoFile) then ("--info-file= '" + infoFile + "'") else ""} \
         ~{if defined(wildcardFile) then ("--wildcard-file= '" + wildcardFile + "'") else ""} \
         ~{if defined(tooShortOutput) then ("--too-short-output= '" + tooShortOutput + "'") else ""} \
         ~{if defined(tooLongOutput) then ("--too-long-output= '" + tooLongOutput + "'") else ""} \
         ~{if defined(untrimmedOutput) then ("--untrimmed-output= '" + untrimmedOutput + "'") else ""} \
         ~{if defined(untrimmedPairedOutput) then ("--untrimmed-paired-output= '" + untrimmedPairedOutput + "'") else ""} \
         ~{if defined(tooShortPairedOutput) then ("--too-short-paired-output= '" + tooShortPairedOutput + "'") else ""} \
         ~{if defined(tooLongPairedOutput) then ("--too-long-paired-output= '" + tooLongPairedOutput + "'") else ""} \
         ~{if defined(inputFileFormat) then ("-f '" + inputFileFormat + "'") else ""} \
         ~{if defined(select_first([cores, 0])) then ("-j " + select_first([cores, 0])) else ''} \
         ~{if defined(adapter_g) then ("-g '" + adapter_g + "'") else ""} \
         ~{if defined(adapter_both) then ("-b '" + adapter_both + "'") else ""} \
         ~{if defined(maximumErrorRate) then ("-e " + maximumErrorRate) else ''} \
         ~{if defined(removeNAdapters) then ("-n " + removeNAdapters) else ''} \
         ~{if defined(overlapRequirement) then ("-O " + overlapRequirement) else ''} \
         ~{if defined(removeNBases) then ("-u " + removeNBases) else ''} \
         ~{if defined(qualityCutoff) then ("-q " + qualityCutoff) else ''} \
         ~{if defined(shortenReadsToLength) then ("-l " + shortenReadsToLength) else ''} \
         ~{if defined(readNamesPrefix) then ("-x '" + readNamesPrefix + "'") else ""} \
         ~{if defined(readNamesSuffix) then ("-y '" + readNamesSuffix + "'") else ""} \
         ~{if defined(minReadLength) then ("-m " + minReadLength) else ''} \
         ~{if defined(maxReadsLength) then ("-M " + maxReadsLength) else ''} \
         ~{if defined(middleReadMatchFile) then ("-r '" + middleReadMatchFile + "'") else ""} \
         ~{if defined(removeMiddle3Adapter) then ("-A '" + removeMiddle3Adapter + "'") else ""} \
         ~{if defined(removeMiddle5Adapter) then ("-G '" + removeMiddle5Adapter + "'") else ""} \
         ~{if defined(removeMiddleBothAdapter) then ("-B '" + removeMiddleBothAdapter + "'") else ""} \
         ~{if defined(removeNBasesFromSecondRead) then ("-U " + removeNBasesFromSecondRead) else ''} \
         ~{if (defined(noMatchAdapterWildcards) && select_first([noMatchAdapterWildcards])) then "-N" else ""} \
         ~{if (defined(colorspace) && select_first([colorspace])) then "-c" else ""} \
         ~{if (defined(doubleEncode) && select_first([doubleEncode])) then "-d" else ""} \
         ~{if (defined(trimPrimer) && select_first([trimPrimer])) then "-t" else ""} \
         ~{if (defined(zeroCap) && select_first([zeroCap])) then "-z" else ""} \
         ~{if length(fastq) > 0 then "'" + sep("' '", fastq) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 5, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/cutadapt:1.18--py37h14c3975_1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4, 4])}G"
       preemptible: 2
     }
     output {
       Array[File] out = glob("*.fastq.gz")
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Cutadapt
   doc: |2-

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

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/cutadapt:1.18--py37h14c3975_1

   inputs:
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
     - string
     - 'null'
     inputBinding:
       prefix: -a
   - id: outputFilename
     label: outputFilename
     doc: |-
       Write trimmed reads to FILE. FASTQ or FASTA format is chosen depending on input. The summary report is sent to standard output. Use '{name}' in FILE to demultiplex reads into multiple files. Default: write to standard output
     type:
     - string
     - 'null'
     default: generated--R1.fastq.gz
     inputBinding:
       prefix: -o
   - id: secondReadFile
     label: secondReadFile
     doc: Write second read in a pair to FILE.
     type:
     - string
     - 'null'
     default: generated--R2.fastq.gz
     inputBinding:
       prefix: -p
   - id: debug
     label: debug
     doc: Print debugging information.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --debug
   - id: noIndels
     label: noIndels
     doc: 'Allow only mismatches in alignments. Default: allow both mismatches and indels'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-indels
   - id: matchReadWildcards
     label: matchReadWildcards
     doc: 'Interpret IUPAC wildcards in reads. Default: False'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --match-read-wildcards
   - id: trimN
     label: trimN
     doc: Trim N's on ends of reads.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --trim-n
   - id: discardCasava
     label: discardCasava
     doc: Discard reads that did not pass CASAVA filtering (header has :Y:).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --discard-casava
   - id: quiet
     label: quiet
     doc: Print only error messages.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --quiet
   - id: stripF3
     label: stripF3
     doc: Strip the _F3 suffix of read names
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --strip-f3
   - id: noZeroCap
     label: noZeroCap
     doc: Disable zero capping
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-zero-cap
   - id: interleaved
     label: interleaved
     doc: Read and write interleaved paired-end reads.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --interleaved
   - id: discardTrimmed
     label: discardTrimmed
     doc: |-
       Discard reads that contain an adapter. Also use -O to avoid discarding too many randomly matching reads!
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --discard-trimmed
   - id: discardUntrimmed
     label: discardUntrimmed
     doc: Discard reads that do not contain an adapter.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --discard-untrimmed
   - id: maq
     label: maq
     doc: |-
       MAQ- and BWA-compatible colorspace output. This enables -c, -d, -t, --strip-f3 and -y '/1'.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --maq
   - id: pairFilter
     label: pairFilter
     doc: |-
       (any|both|first) Which of the reads in a paired-end read have to match the filtering criterion in order for the pair to be filtered. Default: any
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --pair-filter=
   - id: nextseqTrim
     label: nextseqTrim
     doc: |-
       NextSeq-specific quality trimming (each read). Trims also dark cycles appearing as high-quality G bases.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --nextseq-trim=
   - id: action
     label: action
     doc: |-
       What to do with found adapters. trim: remove; mask: replace with 'N' characters; none: leave unchanged (useful with --discard-untrimmed). Default: trim
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --action=
   - id: qualityBase
     label: qualityBase
     doc: |-
       Assume that quality values in FASTQ are encoded as ascii(quality + N). This needs to be set to 64 for some old Illumina FASTQ files. Default: 33
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --quality-base=
   - id: lengthTag
     label: lengthTag
     doc: |-
       Search for TAG followed by a decimal number in the description field of the read. Replace the decimal number with the correct length of the trimmed read. For example, use --length-tag 'length=' to correct fields like 'length=123'.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --length-tag=
   - id: stripSuffix
     label: stripSuffix
     doc: Remove this suffix from read names if present. Can be given multiple times.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --strip-suffix=
   - id: maxN
     label: maxN
     doc: |-
       Discard reads with more than COUNT 'N' bases. If COUNT is a number between 0 and 1, it is interpreted as a fraction of the read length.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-n=
   - id: report
     label: report
     doc: 'Which type of report to print. Default: full'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --report=
   - id: infoFile
     label: infoFile
     doc: |-
       Write information about each read and its adapter matches into FILE. See the documentation for the file format.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --info-file=
   - id: wildcardFile
     label: wildcardFile
     doc: |-
       When the adapter has N wildcard bases, write adapter bases matching wildcard positions to FILE. (Inaccurate with indels.)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --wildcard-file=
   - id: tooShortOutput
     label: tooShortOutput
     doc: |-
       Write reads that are too short (according to length specified by -m) to FILE. Default: discard reads
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --too-short-output=
   - id: tooLongOutput
     label: tooLongOutput
     doc: |-
       Write reads that are too long (according to length specified by -M) to FILE. Default: discard reads
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --too-long-output=
   - id: untrimmedOutput
     label: untrimmedOutput
     doc: |-
       Write reads that do not contain any adapter to FILE. Default: output to same file as trimmed reads
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --untrimmed-output=
   - id: untrimmedPairedOutput
     label: untrimmedPairedOutput
     doc: |-
       Write second read in a pair to this FILE when no adapter was found. Use with --untrimmed-output. Default: output to same file as trimmed reads
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --untrimmed-paired-output=
   - id: tooShortPairedOutput
     label: tooShortPairedOutput
     doc: |-
       Write second read in a pair to this file if pair is too short. Use also --too-short-output.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --too-short-paired-output=
   - id: tooLongPairedOutput
     label: tooLongPairedOutput
     doc: |-
       Write second read in a pair to this file if pair is too long. Use also --too-long-output.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --too-long-paired-output=
   - id: inputFileFormat
     label: inputFileFormat
     doc: |-
       Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. Ignored when reading csfasta/qual files. Default: auto-detect from file name extension.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -f
   - id: cores
     label: cores
     doc: 'Number of CPU cores to use. Use 0 to auto-detect. Default: 1'
     type: int
     default: 0
     inputBinding:
       prefix: -j
   - id: adapter_g
     label: adapter_g
     doc: |-
       Sequence of an adapter ligated to the 5' end (paired data: of the first read). The adapter and any preceding bases are trimmed. Partial matches at the 5' end are allowed. If a '^' character is prepended ('anchoring'), the adapter is only found if it is a prefix of the read.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -g
   - id: adapter_both
     label: adapter_both
     doc: |-
       Sequence of an adapter that may be ligated to the 5' or 3' end (paired data: of the first read). Both types of matches as described under -a und -g are allowed. If the first base of the read is part of the match, the behavior is as with -g, otherwise as with -a. This option is mostly for rescuing failed library preparations - do not use if you know which end your adapter was ligated to!
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -b
   - id: maximumErrorRate
     label: maximumErrorRate
     doc: |-
       Maximum allowed error rate as value between 0 and 1 (no. of errors divided by length of matching region). Default: 0.1 (=10%)
     type:
     - float
     - 'null'
     inputBinding:
       prefix: -e
   - id: removeNAdapters
     label: removeNAdapters
     doc: 'Remove up to COUNT adapters from each read. Default: 1'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -n
   - id: overlapRequirement
     label: overlapRequirement
     doc: |-
       Require MINLENGTH overlap between read and adapter for an adapter to be found. Default: 3
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -O
   - id: removeNBases
     label: removeNBases
     doc: |-
       Remove bases from each read (first read only if paired). If LENGTH is positive, remove bases from the beginning. If LENGTH is negative, remove bases from the end. Can be used twice if LENGTHs have different signs. This is applied *before* adapter trimming.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -u
   - id: qualityCutoff
     label: qualityCutoff
     doc: |-
       --quality-cutoff=[5'CUTOFF,]3'CUTOFF Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. Applied to both reads if data is paired. If one value is given, only the 3' end is trimmed. If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -q
   - id: shortenReadsToLength
     label: shortenReadsToLength
     doc: |-
       Shorten reads to LENGTH. Positive values remove bases at the end while negative ones remove bases at the beginning. This and the following modifications are applied after adapter trimming.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -l
   - id: readNamesPrefix
     label: readNamesPrefix
     doc: |-
       Add this prefix to read names. Use {name} to insert the name of the matching adapter.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -x
   - id: readNamesSuffix
     label: readNamesSuffix
     doc: Add this suffix to read names; can also include {name}
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -y
   - id: minReadLength
     label: minReadLength
     doc: '--minimum-length=LEN[:LEN2] Discard reads shorter than LEN. Default: 0'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -m
   - id: maxReadsLength
     label: maxReadsLength
     doc: '--maximum-length=LEN[:LEN2] Discard reads longer than LEN. Default: no limit'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -M
   - id: middleReadMatchFile
     label: middleReadMatchFile
     doc: |-
       When the adapter matches in the middle of a read, write the rest (after the adapter) to FILE.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -r
   - id: removeMiddle3Adapter
     label: removeMiddle3Adapter
     doc: 3' adapter to be removed from second read in a pair.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -A
   - id: removeMiddle5Adapter
     label: removeMiddle5Adapter
     doc: 5' adapter to be removed from second read in a pair.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -G
   - id: removeMiddleBothAdapter
     label: removeMiddleBothAdapter
     doc: 5'/3 adapter to be removed from second read in a pair.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -B
   - id: removeNBasesFromSecondRead
     label: removeNBasesFromSecondRead
     doc: Remove LENGTH bases from second read in a pair.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -U
   - id: noMatchAdapterWildcards
     label: noMatchAdapterWildcards
     doc: Do not interpret IUPAC wildcards in adapters.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -N
   - id: colorspace
     label: colorspace
     doc: Enable colorspace mode
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -c
   - id: doubleEncode
     label: doubleEncode
     doc: Double-encode colors (map 0,1,2,3,4 to A,C,G,T,N).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -d
   - id: trimPrimer
     label: trimPrimer
     doc: Trim primer base and the first color
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -t
   - id: zeroCap
     label: zeroCap
     doc: |-
       Change negative quality values to zero. Enabled by default in colorspace mode since many tools have problems with negative qualities
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -z

   outputs:
   - id: out
     label: out
     type:
       type: array
       items: File
     outputBinding:
       glob: '*.fastq.gz'
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


