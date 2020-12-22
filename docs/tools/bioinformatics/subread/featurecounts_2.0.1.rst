:orphan:

featureCounts
=============

``featureCounts`` · *1 contributor · 1 version*

FeatureCounts: A General-Purpose Read Summarization Function
This function assigns mapped sequencing reads to genomic features


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.subread.featurecounts.versions import FeatureCounts_2_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "featurecounts_step",
           FeatureCounts_2_0_1(
               bam=None,
               annotationFile=None,
           )
       )
       wf.output("out", source=featurecounts_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for featureCounts:

.. code-block:: bash

   # user inputs
   janis inputs featureCounts > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       annotationFile: annotationFile
       bam:
       - bam_0.bam
       - bam_1.bam




5. Run featureCounts with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       featureCounts





Information
------------

:ID: ``featureCounts``
:URL: `https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts <https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts>`_
:Versions: 2.0.1
:Container: quay.io/biocontainers/subread:2.0.1--hed695b0_0
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-07-16
:Updated: 2020-07-16


Outputs
-----------

======  ========  ===============
name    type      documentation
======  ========  ===============
out     TextFile
======  ========  ===============


Additional configuration (inputs)
---------------------------------

==================  =======================  ====================  ==========  ==================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                type                     prefix                  position  documentation
==================  =======================  ====================  ==========  ==================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
bam                 Array<BAM>                                             10  A list of SAM or BAM format files. They can be either name or location sorted. If no files provided, <stdin> input is expected. Location-sorted paired-end reads are automatically sorted by read names.
annotationFile      File                     -a                                Name of an annotation file. GTF/GFF format by default. See -F option for more format information. Inbuilt annotations (SAF format) is available in 'annotation' directory of the package. Gzipped file is also accepted.
format              Optional<String>         -F                                Specify format of the provided annotation file. Acceptable formats include 'GTF' (or compatible GFF format) and 'SAF'. 'GTF' by default.  For SAF format, please refer to Users Guide.
featureType         Optional<Array<String>>  -t                                Specify feature type(s) in a GTF annotation. If multiple types are provided, they should be separated by ',' with no space in between. 'exon' by default. Rows in the annotation with a matched feature will be extracted and used for read mapping.
attributeType       Optional<String>         -g                                Specify attribute type in GTF annotation. 'gene_id' by default. Meta-features used for read counting will be extracted from annotation using the provided value.
extraAttributes     Optional<Array<String>>  --extraAttributes                 Extract extra attribute types from the provided GTF annotation and include them in the counting output. These attribute types will not be used to group features. If more than one attribute type is provided they should be separated by comma.
chromsomeAlias      Optional<String>         -A                                Provide a chromosome name alias file to match chr names inannotation with those in the reads. This should be a two-column comma-delimited text file. Its first column should include chr names in the annotation and its second column should include chr names in the reads. Chr names are case sensitive. No column header should be included in the file.
featureLevel        Optional<Boolean>        -f                                Perform read counting at feature level (eg. counting reads for exons rather than genes).
overlap             Optional<Boolean>        -O                                Assign reads to all their overlapping meta-features (or features if -f is specified).
minOverlap          Optional<Integer>        --minOverlap                      Minimum number of overlapping bases in a read that isrequired for read assignment. 1 by default. Number ofoverlapping bases is counted from both reads if pairedend. If a negative value is provided, then a gap of upto specified size will be allowed between read and the feature that the read is assigned to.
fracOverlap         Optional<Float>          --fracOverlap                     Minimum fraction of overlapping bases in a read that isrequired for read assignment. Value should be within range [0,1]. 0 by default. Number of overlapping bases is counted from both reads if paired end. Both this option and '--minOverlap' option need to be satisfied for read assignment.
fracOverlapFeature  Optional<Float>          --fracOverlapFeature              Minimum fraction of overlapping bases in a feature that is required for read assignment. Value should be within range [0,1]. 0 by default.
largestOverlap      Optional<Boolean>        --largestOverlap                  Assign reads to a meta-feature/feature that has the  largest number of overlapping bases.
nonOverlap          Optional<Integer>        --nonOverlap                      Maximum number of non-overlapping bases in a read (or a read pair) that is allowed when being assigned to a feature. No limit is set by default.
nonOverlapFeature   Optional<Integer>        --nonOverlapFeature               Maximum number of non-overlapping bases in a feature that is allowed in read assignment. No limit is set by default.
readExtensionFive   Optional<Integer>        --readExtension5                  Reads are extended upstream by <int> bases from their 5' end.
readExtensionThree  Optional<String>         --readExtension3                  Reads are extended upstream by <int> bases from their 3' end.
readToPos           Optional<String>         --read2pos                        Reduce reads to their 5' most base or 3' most base. Read counting is then performed based on the single base the read is reduced to.
multiMapping        Optional<Boolean>        -M                                Multi-mapping reads will also be counted. For a multi-mapping read, all its reported alignments will be counted. The 'NH' tag in BAM/SAM input is used to detect multi-mapping reads.
fration             Optional<Boolean>        --fraction                        Assign fractional counts to features. This option must be used together with '-M' or '-O' or both. When '-M' is specified, each reported alignment from a multi-mapping read (identified via 'NH' tag) will carry a fractional count of 1/x, instead of 1 (one), where x is the total number of alignments reported for the same read. When '-O' is specified, each overlapping feature will receive a fractional count of 1/y, where y is the total number of features overlapping with the read. When both '-M' and '-O' are specified, each alignment will carry a fractional count of 1/(x*y).
quality             Optional<String>         -Q                                The minimum mapping quality score a read must satisfy in order to be counted. For paired-end reads, at least one end should satisfy this criteria. 0 by default.
splitOnly           Optional<Boolean>        --splitOnly                       Count split alignments only (ie. alignments with CIGAR string containing 'N'). An example of split alignments is exon-spanning reads in RNA-seq data.
nonSplitOnly        Optional<Boolean>        --nonSplitOnly                    If specified, only non-split alignments (CIGAR strings do not contain letter 'N') will be counted. All the other alignments will be ignored.
primary             Optional<Boolean>        --primary                         Count primary alignments only. Primary alignments are identified using bit 0x100 in SAM/BAM FLAG field.
ignoreDup           Optional<Boolean>        --ignoreDup                       Ignore duplicate reads in read counting. Duplicate reads are identified using bit Ox400 in BAM/SAM FLAG field. The whole read pair is ignored if one of the reads is a duplicate read for paired end data.
strandness          Optional<String>         -                                 Perform strand-specific read counting. A single integer value (applied to all input files) or a string of comma-separated values (applied to each corresponding input file) should be provided. Possible values include: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). Default value is 0 (ie. unstranded read counting carried out for all input files).
junction            Optional<String>         -J                                Count number of reads supporting each exon-exon junction. Junctions were identified from those exon-spanning reads in the input (containing 'N' in CIGAR string). Counting results are saved to a file named '<output_file>.jcounts'
genome              Optional<File>           -G                                Provide the name of a FASTA-format file that contains thereference sequences used in read mapping that produced the provided SAM/BAM files. This optional argument can be used with '-J' option to improve read counting for junctions.
pairEnd             Optional<Boolean>        -p                                If specified, fragments (or templates) will be counted instead of reads. This option is only applicable for paired-end reads; single-end reads are always counted as reads.
both                Optional<Boolean>        -B                                Only count read pairs that have both ends aligned.
pairEndDistance     Optional<Boolean>        -P                                Check validity of paired-end distance when counting read  pairs. Use -d and -D to set thresholds.
minDistance         Optional<Integer>        -d                                Minimum fragment/template length, 50 by default.
maxDistance         Optional<Integer>        -D                                Maximum fragment/template length, 600 by default.
countRead           Optional<Boolean>        -C                                Do not count read pairs that have their two ends mapping to different chromosomes or mapping to same chromosome but on different strands.
doNotSort           Optional<Boolean>        --donotsort                       Do not sort reads in BAM/SAM input. Note that reads from the same pair are required to be located next to each other in the input.
threads             Optional<Integer>        -T                                Number of the threads. 1 by default.
byReadGroup         Optional<Boolean>        --byReadGroup                     Assign reads by read group. 'RG' tag is required to be present in the input BAM/SAM files.
longRead            Optional<Boolean>        -L                                Count long reads such as Nanopore and PacBio reads. Long read counting can only run in one thread and only reads (not read-pairs) can be counted. There is no limitation on the number of 'M' operations allowed in a CIGAR string in long read counting.
outputFormat        Optional<String>         -R                                Output detailed assignment results for each read or read-pair. Results are saved to a file that is in one of the following formats: CORE, SAM and BAM. See Users Guide for more info about these formats.
outputDirectory     Optional<String>         --Rpath                           Specify a directory to save the detailed assignment results. If unspecified, the directory where counting results are saved is used.
tmpDir              Optional<String>         --tmpDir                          Directory under which intermediate files are saved (later removed). By default, intermediate files will be saved to the directory specified in '-o' argument.
maxMOp              Optional<Integer>        --maxMOp                          Maximum number of 'M' operations allowed in a CIGAR string. 10 by default. Both 'X' and '=' are treated as 'M' and adjacent 'M' operations are merged in the CIGAR string.
outputFilename      Optional<Filename>       -o                                Name of output file including read counts. A separate file including summary statistics of counting results is also included in the output ('<string>.summary'). Both files are in tab delimited format.
==================  =======================  ====================  ==========  ==================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task featureCounts {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       String? format
       Array[String]? featureType
       String? attributeType
       Array[String]? extraAttributes
       String? chromsomeAlias
       Boolean? featureLevel
       Boolean? overlap
       Int? minOverlap
       Float? fracOverlap
       Float? fracOverlapFeature
       Boolean? largestOverlap
       Int? nonOverlap
       Int? nonOverlapFeature
       Int? readExtensionFive
       String? readExtensionThree
       String? readToPos
       Boolean? multiMapping
       Boolean? fration
       String? quality
       Boolean? splitOnly
       Boolean? nonSplitOnly
       Boolean? primary
       Boolean? ignoreDup
       String? strandness
       String? junction
       File? genome
       Boolean? pairEnd
       Boolean? both
       Boolean? pairEndDistance
       Int? minDistance
       Int? maxDistance
       Boolean? countRead
       Boolean? doNotSort
       Int? threads
       Boolean? byReadGroup
       Boolean? longRead
       String? outputFormat
       String? outputDirectory
       String? tmpDir
       Int? maxMOp
       Array[File] bam
       String? outputFilename
       File annotationFile
     }
     command <<<
       set -e
        featureCounts \
         ~{if defined(format) then ("-F '" + format + "'") else ""} \
         ~{if (defined(featureType) && length(select_first([featureType])) > 0) then "-t '" + sep("','", select_first([featureType])) + "'" else ""} \
         ~{if defined(attributeType) then ("-g '" + attributeType + "'") else ""} \
         ~{if (defined(extraAttributes) && length(select_first([extraAttributes])) > 0) then "--extraAttributes '" + sep("','", select_first([extraAttributes])) + "'" else ""} \
         ~{if defined(chromsomeAlias) then ("-A '" + chromsomeAlias + "'") else ""} \
         ~{if (defined(featureLevel) && select_first([featureLevel])) then "-f" else ""} \
         ~{if (defined(overlap) && select_first([overlap])) then "-O" else ""} \
         ~{if defined(minOverlap) then ("--minOverlap " + minOverlap) else ''} \
         ~{if defined(fracOverlap) then ("--fracOverlap " + fracOverlap) else ''} \
         ~{if defined(fracOverlapFeature) then ("--fracOverlapFeature " + fracOverlapFeature) else ''} \
         ~{if (defined(largestOverlap) && select_first([largestOverlap])) then "--largestOverlap" else ""} \
         ~{if defined(nonOverlap) then ("--nonOverlap " + nonOverlap) else ''} \
         ~{if defined(nonOverlapFeature) then ("--nonOverlapFeature " + nonOverlapFeature) else ''} \
         ~{if defined(readExtensionFive) then ("--readExtension5 " + readExtensionFive) else ''} \
         ~{if defined(readExtensionThree) then ("--readExtension3 '" + readExtensionThree + "'") else ""} \
         ~{if defined(readToPos) then ("--read2pos '" + readToPos + "'") else ""} \
         ~{if (defined(multiMapping) && select_first([multiMapping])) then "-M" else ""} \
         ~{if (defined(fration) && select_first([fration])) then "--fraction" else ""} \
         ~{if defined(quality) then ("-Q '" + quality + "'") else ""} \
         ~{if (defined(splitOnly) && select_first([splitOnly])) then "--splitOnly" else ""} \
         ~{if (defined(nonSplitOnly) && select_first([nonSplitOnly])) then "--nonSplitOnly" else ""} \
         ~{if (defined(primary) && select_first([primary])) then "--primary" else ""} \
         ~{if (defined(ignoreDup) && select_first([ignoreDup])) then "--ignoreDup" else ""} \
         ~{if defined(strandness) then ("- '" + strandness + "'") else ""} \
         ~{if defined(junction) then ("-J '" + junction + "'") else ""} \
         ~{if defined(genome) then ("-G '" + genome + "'") else ""} \
         ~{if (defined(pairEnd) && select_first([pairEnd])) then "-p" else ""} \
         ~{if (defined(both) && select_first([both])) then "-B" else ""} \
         ~{if (defined(pairEndDistance) && select_first([pairEndDistance])) then "-P" else ""} \
         ~{if defined(minDistance) then ("-d " + minDistance) else ''} \
         ~{if defined(maxDistance) then ("-D " + maxDistance) else ''} \
         ~{if (defined(countRead) && select_first([countRead])) then "-C" else ""} \
         ~{if (defined(doNotSort) && select_first([doNotSort])) then "--donotsort" else ""} \
         ~{if defined(threads) then ("-T " + threads) else ''} \
         ~{if (defined(byReadGroup) && select_first([byReadGroup])) then "--byReadGroup" else ""} \
         ~{if (defined(longRead) && select_first([longRead])) then "-L" else ""} \
         ~{if defined(outputFormat) then ("-R '" + outputFormat + "'") else ""} \
         ~{if defined(outputDirectory) then ("--Rpath '" + outputDirectory + "'") else ""} \
         ~{if defined(tmpDir) then ("--tmpDir '" + tmpDir + "'") else ""} \
         ~{if defined(maxMOp) then ("--maxMOp " + maxMOp) else ''} \
         -o '~{select_first([outputFilename, "generated.txt"])}' \
         -a '~{annotationFile}' \
         ~{if length(bam) > 0 then "'" + sep("' '", bam) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/subread:2.0.1--hed695b0_0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.txt"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: featureCounts
   doc: |-
     FeatureCounts: A General-Purpose Read Summarization Function
     This function assigns mapped sequencing reads to genomic features

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/subread:2.0.1--hed695b0_0

   inputs:
   - id: format
     label: format
     doc: |-
       Specify format of the provided annotation file. Acceptable formats include 'GTF' (or compatible GFF format) and 'SAF'. 'GTF' by default.  For SAF format, please refer to Users Guide.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -F
   - id: featureType
     label: featureType
     doc: |-
       Specify feature type(s) in a GTF annotation. If multiple types are provided, they should be separated by ',' with no space in between. 'exon' by default. Rows in the annotation with a matched feature will be extracted and used for read mapping.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: -t
       itemSeparator: ','
   - id: attributeType
     label: attributeType
     doc: |-
       Specify attribute type in GTF annotation. 'gene_id' by default. Meta-features used for read counting will be extracted from annotation using the provided value.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -g
   - id: extraAttributes
     label: extraAttributes
     doc: |-
       Extract extra attribute types from the provided GTF annotation and include them in the counting output. These attribute types will not be used to group features. If more than one attribute type is provided they should be separated by comma.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --extraAttributes
       itemSeparator: ','
   - id: chromsomeAlias
     label: chromsomeAlias
     doc: |-
       Provide a chromosome name alias file to match chr names inannotation with those in the reads. This should be a two-column comma-delimited text file. Its first column should include chr names in the annotation and its second column should include chr names in the reads. Chr names are case sensitive. No column header should be included in the file.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -A
   - id: featureLevel
     label: featureLevel
     doc: |-
       Perform read counting at feature level (eg. counting reads for exons rather than genes).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -f
   - id: overlap
     label: overlap
     doc: |-
       Assign reads to all their overlapping meta-features (or features if -f is specified).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -O
   - id: minOverlap
     label: minOverlap
     doc: |-
       Minimum number of overlapping bases in a read that isrequired for read assignment. 1 by default. Number ofoverlapping bases is counted from both reads if pairedend. If a negative value is provided, then a gap of upto specified size will be allowed between read and the feature that the read is assigned to.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --minOverlap
   - id: fracOverlap
     label: fracOverlap
     doc: |-
       Minimum fraction of overlapping bases in a read that isrequired for read assignment. Value should be within range [0,1]. 0 by default. Number of overlapping bases is counted from both reads if paired end. Both this option and '--minOverlap' option need to be satisfied for read assignment.
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --fracOverlap
   - id: fracOverlapFeature
     label: fracOverlapFeature
     doc: |-
       Minimum fraction of overlapping bases in a feature that is required for read assignment. Value should be within range [0,1]. 0 by default.
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --fracOverlapFeature
   - id: largestOverlap
     label: largestOverlap
     doc: |-
       Assign reads to a meta-feature/feature that has the  largest number of overlapping bases.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --largestOverlap
   - id: nonOverlap
     label: nonOverlap
     doc: |-
       Maximum number of non-overlapping bases in a read (or a read pair) that is allowed when being assigned to a feature. No limit is set by default.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --nonOverlap
   - id: nonOverlapFeature
     label: nonOverlapFeature
     doc: |-
       Maximum number of non-overlapping bases in a feature that is allowed in read assignment. No limit is set by default.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --nonOverlapFeature
   - id: readExtensionFive
     label: readExtensionFive
     doc: Reads are extended upstream by <int> bases from their 5' end.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --readExtension5
   - id: readExtensionThree
     label: readExtensionThree
     doc: Reads are extended upstream by <int> bases from their 3' end.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --readExtension3
   - id: readToPos
     label: readToPos
     doc: |-
       Reduce reads to their 5' most base or 3' most base. Read counting is then performed based on the single base the read is reduced to.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read2pos
   - id: multiMapping
     label: multiMapping
     doc: |-
       Multi-mapping reads will also be counted. For a multi-mapping read, all its reported alignments will be counted. The 'NH' tag in BAM/SAM input is used to detect multi-mapping reads.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -M
   - id: fration
     label: fration
     doc: |-
       Assign fractional counts to features. This option must be used together with '-M' or '-O' or both. When '-M' is specified, each reported alignment from a multi-mapping read (identified via 'NH' tag) will carry a fractional count of 1/x, instead of 1 (one), where x is the total number of alignments reported for the same read. When '-O' is specified, each overlapping feature will receive a fractional count of 1/y, where y is the total number of features overlapping with the read. When both '-M' and '-O' are specified, each alignment will carry a fractional count of 1/(x*y).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --fraction
   - id: quality
     label: quality
     doc: |-
       The minimum mapping quality score a read must satisfy in order to be counted. For paired-end reads, at least one end should satisfy this criteria. 0 by default.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -Q
   - id: splitOnly
     label: splitOnly
     doc: |-
       Count split alignments only (ie. alignments with CIGAR string containing 'N'). An example of split alignments is exon-spanning reads in RNA-seq data.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --splitOnly
   - id: nonSplitOnly
     label: nonSplitOnly
     doc: |-
       If specified, only non-split alignments (CIGAR strings do not contain letter 'N') will be counted. All the other alignments will be ignored.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --nonSplitOnly
   - id: primary
     label: primary
     doc: |-
       Count primary alignments only. Primary alignments are identified using bit 0x100 in SAM/BAM FLAG field.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --primary
   - id: ignoreDup
     label: ignoreDup
     doc: |-
       Ignore duplicate reads in read counting. Duplicate reads are identified using bit Ox400 in BAM/SAM FLAG field. The whole read pair is ignored if one of the reads is a duplicate read for paired end data.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ignoreDup
   - id: strandness
     label: strandness
     doc: |-
       Perform strand-specific read counting. A single integer value (applied to all input files) or a string of comma-separated values (applied to each corresponding input file) should be provided. Possible values include: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). Default value is 0 (ie. unstranded read counting carried out for all input files).
     type:
     - string
     - 'null'
     inputBinding:
       prefix: '-'
   - id: junction
     label: junction
     doc: |-
       Count number of reads supporting each exon-exon junction. Junctions were identified from those exon-spanning reads in the input (containing 'N' in CIGAR string). Counting results are saved to a file named '<output_file>.jcounts'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -J
   - id: genome
     label: genome
     doc: |-
       Provide the name of a FASTA-format file that contains thereference sequences used in read mapping that produced the provided SAM/BAM files. This optional argument can be used with '-J' option to improve read counting for junctions.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -G
   - id: pairEnd
     label: pairEnd
     doc: |-
       If specified, fragments (or templates) will be counted instead of reads. This option is only applicable for paired-end reads; single-end reads are always counted as reads.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -p
   - id: both
     label: both
     doc: Only count read pairs that have both ends aligned.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -B
   - id: pairEndDistance
     label: pairEndDistance
     doc: |-
       Check validity of paired-end distance when counting read  pairs. Use -d and -D to set thresholds.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -P
   - id: minDistance
     label: minDistance
     doc: Minimum fragment/template length, 50 by default.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -d
   - id: maxDistance
     label: maxDistance
     doc: Maximum fragment/template length, 600 by default.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -D
   - id: countRead
     label: countRead
     doc: |-
       Do not count read pairs that have their two ends mapping to different chromosomes or mapping to same chromosome but on different strands.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -C
   - id: doNotSort
     label: doNotSort
     doc: |-
       Do not sort reads in BAM/SAM input. Note that reads from the same pair are required to be located next to each other in the input.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --donotsort
   - id: threads
     label: threads
     doc: Number of the threads. 1 by default.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -T
   - id: byReadGroup
     label: byReadGroup
     doc: |-
       Assign reads by read group. 'RG' tag is required to be present in the input BAM/SAM files.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --byReadGroup
   - id: longRead
     label: longRead
     doc: |-
       Count long reads such as Nanopore and PacBio reads. Long read counting can only run in one thread and only reads (not read-pairs) can be counted. There is no limitation on the number of 'M' operations allowed in a CIGAR string in long read counting.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -L
   - id: outputFormat
     label: outputFormat
     doc: |-
       Output detailed assignment results for each read or read-pair. Results are saved to a file that is in one of the following formats: CORE, SAM and BAM. See Users Guide for more info about these formats.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -R
   - id: outputDirectory
     label: outputDirectory
     doc: |-
       Specify a directory to save the detailed assignment results. If unspecified, the directory where counting results are saved is used.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --Rpath
   - id: tmpDir
     label: tmpDir
     doc: |-
       Directory under which intermediate files are saved (later removed). By default, intermediate files will be saved to the directory specified in '-o' argument.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --tmpDir
   - id: maxMOp
     label: maxMOp
     doc: |-
       Maximum number of 'M' operations allowed in a CIGAR string. 10 by default. Both 'X' and '=' are treated as 'M' and adjacent 'M' operations are merged in the CIGAR string.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --maxMOp
   - id: bam
     label: bam
     doc: |-
       A list of SAM or BAM format files. They can be either name or location sorted. If no files provided, <stdin> input is expected. Location-sorted paired-end reads are automatically sorted by read names.
     type:
       type: array
       items: File
     inputBinding:
       position: 10
   - id: outputFilename
     label: outputFilename
     doc: |-
       Name of output file including read counts. A separate file including summary statistics of counting results is also included in the output ('<string>.summary'). Both files are in tab delimited format.
     type:
     - string
     - 'null'
     default: generated.txt
     inputBinding:
       prefix: -o
   - id: annotationFile
     label: annotationFile
     doc: |-
       Name of an annotation file. GTF/GFF format by default. See -F option for more format information. Inbuilt annotations (SAF format) is available in 'annotation' directory of the package. Gzipped file is also accepted.
     type: File
     inputBinding:
       prefix: -a

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.txt
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - ''
   - featureCounts
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: featureCounts


