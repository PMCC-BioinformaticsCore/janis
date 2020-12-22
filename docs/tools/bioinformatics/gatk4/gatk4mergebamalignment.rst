:orphan:

GATK4: Merge SAM or BAM with unmapped BAM file
=======================================================================

``Gatk4MergeBamAlignment`` · *2 contributors · 4 versions*

Merges SAM/BAM file with an unmapped BAM file


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.mergebamalignment.versions import Gatk4MergeBamAlignment_4_1_4

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4mergebamalignment_step",
           Gatk4MergeBamAlignment_4_1_4(
               ubam=None,
               bam=None,
           )
       )
       wf.output("out", source=gatk4mergebamalignment_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4MergeBamAlignment:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4MergeBamAlignment > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam:
       - bam_0.sam
       - bam_1.sam
       ubam: ubam.bam




5. Run Gatk4MergeBamAlignment with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4MergeBamAlignment





Information
------------

:ID: ``Gatk4MergeBamAlignment``
:URL: `https://gatk.broadinstitute.org/hc/en-us/articles/360037225832-MergeBamAlignment-Picard- <https://gatk.broadinstitute.org/hc/en-us/articles/360037225832-MergeBamAlignment-Picard->`_
:Versions: 4.1.4.1, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.4.1
:Authors: Michael Franklin (@illusional), Matthias De Smet(@matthdsm)
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2020-02-26


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     BAM
======  ======  ===============


Additional configuration (inputs)
---------------------------------

=============================  ==========================  ==================================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
name                           type                        prefix                                position  documentation
=============================  ==========================  ==================================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
ubam                           BAM                         --UNMAPPED_BAM                              10  Original SAM or BAM file of unmapped reads, which must be in queryname order.
bam                            Array<SAM>                  --ALIGNED_BAM                               10  SAM or BAM file(s) with alignment data.
javaOptions                    Optional<Array<String>>
compression_level              Optional<Integer>                                                           Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
reference                      Optional<FastaWithIndexes>  --REFERENCE_SEQUENCE                        10  Reference sequence file.
outputFilename                 Optional<Filename>          --OUTPUT                                    10  Merged SAM or BAM file to write to.
addMateCigar                   Optional<Boolean>           --ADD_MATE_CIGAR                            11  Adds the mate CIGAR tag (MC)
alignedReadsOnly               Optional<Boolean>           --ALIGNED_READS_ONLY                        11  Whether to output only aligned reads.
alignerProperPairFlags         Optional<Boolean>           --ALIGNER_PROPER_PAIR_FLAGS                 11  Use the aligner's idea of what a proper pair is rather than computing in this program.
argumentsFile                  Optional<Array<File>>       --arguments_file                            11  read one or more arguments files and add them to the command line
attributesToRemove             Optional<Array<String>>     --ATTRIBUTES_TO_REMOVE                      11  Attributes from the alignment record that should be removed when merging.
attributesToRetain             Optional<Array<String>>     --ATTRIBUTES_TO_RETAIN                      11  Reserved alignment attributes (tags starting with X, Y, or Z) that should be brought over from the alignment data when merging.
attributesToReverse            Optional<Array<String>>     --ATTRIBUTES_TO_REVERSE                     11  Attributes on negative strand reads that need to be reversed.
attributesToReverseComplement  Optional<Array<String>>     --ATTRIBUTES_TO_REVERSE_COMPLEMENT          11  Attributes on negative strand reads that need to be reverse complemented.
clipAdapter                    Optional<Boolean>           --CLIP_ADAPTERS                             11  Whether to clip adapters where identified.
clipOverlappingReads           Optional<Boolean>           --CLIP_OVERLAPPING_READS                    11  For paired reads, soft clip the 3' end of each read if necessary so that it does not extend past the 5' end of its mate.
expectedOrientations           Optional<Array<String>>     --EXPECTED_ORIENTATIONS                     11  The expected orientation of proper read pairs.
includeSecondaryAlginments     Optional<Boolean>           --INCLUDE_SECONDARY_ALIGNMENTS              11  If false, do not write secondary alignments to output.
isBisulfiteSequencing          Optional<Boolean>           --IS_BISULFITE_SEQUENCE                     11  Whether the lane is bisulfite sequence (used when calculating the NM tag).
matchingDictionaryTags         Optional<Array<String>>     --MATCHING_DICTIONARY_TAGS                  11  List of Sequence Records tags that must be equal (if present) in the reference dictionary and in the aligned file.
maxInsertionsOrDeletions       Optional<Integer>           --MAX_INSERTIONS_OR_DELETIONS               11  The maximum number of insertions or deletions permitted for an alignment to be included.
minUnclippedBases              Optional<Integer>           --MIN_UNCLIPPED_BASES                       11  If UNMAP_CONTAMINANT_READS is set, require this many unclipped bases or else the read will be marked as contaminant.
primaryAlignmentStrategy       Optional<Integer>           --PRIMARY_ALIGNMENT_STRATEGY                11  Strategy for selecting primary alignment when the aligner has provided more than one alignment for a pair or fragment, and none are marked as primary, more than one is marked as primary, or the primary alignment is filtered out for some reason.
programGroupCommandLine        Optional<String>            --PROGRAM_GROUP_COMMAND_LINE                11  The command line of the program group.
programGroupName               Optional<String>            --PROGRAM_GROUP_NAME                        11  The name of the program group.
programGroupVersion            Optional<String>            --PROGRAM_GROUP_VERSION                     11  The version of the program group.
programRecordId                Optional<String>            --PROGRAM_RECORD_ID                         11  The program group ID of the aligner.
sortOrder                      Optional<String>            -SO                                         10  The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
unmapContaminantReads          Optional<Boolean>           --UNMAP_CONTAMINANT_READS                   11  Detect reads originating from foreign organisms (e.g. bacterial DNA in a non-bacterial sample),and unmap + label those reads accordingly.
unmappedReadStrategy           Optional<String>            --UNMAPPED_READ_STRATEGY                    11  How to deal with alignment information in reads that are being unmapped (e.g. due to cross-species contamination.) Currently ignored unless UNMAP_CONTAMINANT_READS = true.
addPgTagToReads                Optional<Boolean>           --ADD_PG_TAG_TO_READS                       11  Add PG tag to each read in a SAM or BAM
compressionLevel               Optional<Integer>           --COMPRESSION_LEVEL                         11  Compression level for all compressed files created (e.g. BAM and GELI).
createIndex                    Optional<Boolean>           --CREATE_INDEX                              11  Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File                  Optional<Boolean>           --CREATE_MD5_FILE                           11  Whether to create an MD5 digest for any BAM or FASTQ files created.
maxRecordsInRam                Optional<Integer>           --MAX_RECORDS_IN_RAM                        11  When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
quiet                          Optional<Boolean>           --QUIET                                     11  Whether to suppress job-summary info on System.err.
tmpDir                         Optional<String>            --TMP_DIR                                   11  Undocumented option
useJdkDeflater                 Optional<Boolean>           --use_jdk_deflater                          11  Whether to use the JdkDeflater (as opposed to IntelDeflater)
useJdkInflater                 Optional<Boolean>           --use_jdk_inflater                          11  Whether to use the JdkInflater (as opposed to IntelInflater)
validationStringency           Optional<String>            --VALIDATION_STRINGENCY                     11  Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
verbosity                      Optional<String>            --verbosity                                 11  The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
=============================  ==========================  ==================================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4MergeBamAlignment {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       File ubam
       Array[File] bam
       File? reference
       File? reference_fai
       File? reference_amb
       File? reference_ann
       File? reference_bwt
       File? reference_pac
       File? reference_sa
       File? reference_dict
       String? outputFilename
       Boolean? addMateCigar
       Boolean? alignedReadsOnly
       Boolean? alignerProperPairFlags
       Array[File]? argumentsFile
       Array[String]? attributesToRemove
       Array[String]? attributesToRetain
       Array[String]? attributesToReverse
       Array[String]? attributesToReverseComplement
       Boolean? clipAdapter
       Boolean? clipOverlappingReads
       Array[String]? expectedOrientations
       Boolean? includeSecondaryAlginments
       Boolean? isBisulfiteSequencing
       Array[String]? matchingDictionaryTags
       Int? maxInsertionsOrDeletions
       Int? minUnclippedBases
       Int? primaryAlignmentStrategy
       String? programGroupCommandLine
       String? programGroupName
       String? programGroupVersion
       String? programRecordId
       String? sortOrder
       Boolean? unmapContaminantReads
       String? unmappedReadStrategy
       Boolean? addPgTagToReads
       Int? compressionLevel
       Boolean? createIndex
       Boolean? createMd5File
       Int? maxRecordsInRam
       Boolean? quiet
       String? tmpDir
       Boolean? useJdkDeflater
       Boolean? useJdkInflater
       String? validationStringency
       String? verbosity
     }
     command <<<
       set -e
       gatk MergeBamAlignment \
         --java-options '-Xmx~{((select_first([runtime_memory, 4, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         --UNMAPPED_BAM '~{ubam}' \
         ~{if length(bam) > 0 then "--ALIGNED_BAM '" + sep("' --ALIGNED_BAM '", bam) + "'" else ""} \
         ~{if defined(reference) then ("--REFERENCE_SEQUENCE '" + reference + "'") else ""} \
         --OUTPUT '~{select_first([outputFilename, "generated.bam"])}' \
         ~{if defined(sortOrder) then ("-SO '" + sortOrder + "'") else ""} \
         ~{if (defined(addMateCigar) && select_first([addMateCigar])) then "--ADD_MATE_CIGAR" else ""} \
         ~{if (defined(alignedReadsOnly) && select_first([alignedReadsOnly])) then "--ALIGNED_READS_ONLY" else ""} \
         ~{if (defined(alignerProperPairFlags) && select_first([alignerProperPairFlags])) then "--ALIGNER_PROPER_PAIR_FLAGS" else ""} \
         ~{if (defined(argumentsFile) && length(select_first([argumentsFile])) > 0) then "--arguments_file '" + sep("' '", select_first([argumentsFile])) + "'" else ""} \
         ~{if (defined(attributesToRemove) && length(select_first([attributesToRemove])) > 0) then "--ATTRIBUTES_TO_REMOVE '" + sep("' '", select_first([attributesToRemove])) + "'" else ""} \
         ~{if (defined(attributesToRetain) && length(select_first([attributesToRetain])) > 0) then "--ATTRIBUTES_TO_RETAIN '" + sep("' '", select_first([attributesToRetain])) + "'" else ""} \
         ~{if (defined(attributesToReverse) && length(select_first([attributesToReverse])) > 0) then "--ATTRIBUTES_TO_REVERSE '" + sep("' '", select_first([attributesToReverse])) + "'" else ""} \
         ~{if (defined(attributesToReverseComplement) && length(select_first([attributesToReverseComplement])) > 0) then "--ATTRIBUTES_TO_REVERSE_COMPLEMENT '" + sep("' '", select_first([attributesToReverseComplement])) + "'" else ""} \
         ~{if (defined(clipAdapter) && select_first([clipAdapter])) then "--CLIP_ADAPTERS" else ""} \
         ~{if (defined(clipOverlappingReads) && select_first([clipOverlappingReads])) then "--CLIP_OVERLAPPING_READS" else ""} \
         ~{if (defined(expectedOrientations) && length(select_first([expectedOrientations])) > 0) then "--EXPECTED_ORIENTATIONS '" + sep("' '", select_first([expectedOrientations])) + "'" else ""} \
         ~{if (defined(includeSecondaryAlginments) && select_first([includeSecondaryAlginments])) then "--INCLUDE_SECONDARY_ALIGNMENTS" else ""} \
         ~{if (defined(isBisulfiteSequencing) && select_first([isBisulfiteSequencing])) then "--IS_BISULFITE_SEQUENCE" else ""} \
         ~{if (defined(matchingDictionaryTags) && length(select_first([matchingDictionaryTags])) > 0) then "--MATCHING_DICTIONARY_TAGS '" + sep("' '", select_first([matchingDictionaryTags])) + "'" else ""} \
         ~{if defined(maxInsertionsOrDeletions) then ("--MAX_INSERTIONS_OR_DELETIONS " + maxInsertionsOrDeletions) else ''} \
         ~{if defined(minUnclippedBases) then ("--MIN_UNCLIPPED_BASES " + minUnclippedBases) else ''} \
         ~{if defined(primaryAlignmentStrategy) then ("--PRIMARY_ALIGNMENT_STRATEGY " + primaryAlignmentStrategy) else ''} \
         ~{if defined(programGroupCommandLine) then ("--PROGRAM_GROUP_COMMAND_LINE '" + programGroupCommandLine + "'") else ""} \
         ~{if defined(programGroupName) then ("--PROGRAM_GROUP_NAME '" + programGroupName + "'") else ""} \
         ~{if defined(programGroupVersion) then ("--PROGRAM_GROUP_VERSION '" + programGroupVersion + "'") else ""} \
         ~{if defined(programRecordId) then ("--PROGRAM_RECORD_ID '" + programRecordId + "'") else ""} \
         ~{if (defined(unmapContaminantReads) && select_first([unmapContaminantReads])) then "--UNMAP_CONTAMINANT_READS" else ""} \
         ~{if defined(unmappedReadStrategy) then ("--UNMAPPED_READ_STRATEGY '" + unmappedReadStrategy + "'") else ""} \
         ~{if (defined(addPgTagToReads) && select_first([addPgTagToReads])) then "--ADD_PG_TAG_TO_READS" else ""} \
         ~{if defined(compressionLevel) then ("--COMPRESSION_LEVEL " + compressionLevel) else ''} \
         ~{if select_first([createIndex, true]) then "--CREATE_INDEX" else ""} \
         ~{if (defined(createMd5File) && select_first([createMd5File])) then "--CREATE_MD5_FILE" else ""} \
         ~{if defined(maxRecordsInRam) then ("--MAX_RECORDS_IN_RAM " + maxRecordsInRam) else ''} \
         ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
         ~{if defined(select_first([tmpDir, "/tmp/"])) then ("--TMP_DIR '" + select_first([tmpDir, "/tmp/"]) + "'") else ""} \
         ~{if (defined(useJdkDeflater) && select_first([useJdkDeflater])) then "--use_jdk_deflater" else ""} \
         ~{if (defined(useJdkInflater) && select_first([useJdkInflater])) then "--use_jdk_inflater" else ""} \
         ~{if defined(validationStringency) then ("--VALIDATION_STRINGENCY '" + validationStringency + "'") else ""} \
         ~{if defined(verbosity) then ("--verbosity '" + verbosity + "'") else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.4.1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.bam"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: Merge SAM or BAM with unmapped BAM file'
   doc: Merges SAM/BAM file with an unmapped BAM file

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.4.1

   inputs:
   - id: javaOptions
     label: javaOptions
     type:
     - type: array
       items: string
     - 'null'
   - id: compression_level
     label: compression_level
     doc: |-
       Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
     type:
     - int
     - 'null'
   - id: ubam
     label: ubam
     doc: Original SAM or BAM file of unmapped reads, which must be in queryname order.
     type: File
     inputBinding:
       prefix: --UNMAPPED_BAM
       position: 10
   - id: bam
     label: bam
     doc: SAM or BAM file(s) with alignment data.
     type:
       type: array
       inputBinding:
         prefix: --ALIGNED_BAM
       items: File
     inputBinding:
       position: 10
   - id: reference
     label: reference
     doc: Reference sequence file.
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .fai
     - pattern: .amb
     - pattern: .ann
     - pattern: .bwt
     - pattern: .pac
     - pattern: .sa
     - pattern: ^.dict
     inputBinding:
       prefix: --REFERENCE_SEQUENCE
       position: 10
   - id: outputFilename
     label: outputFilename
     doc: Merged SAM or BAM file to write to.
     type:
     - string
     - 'null'
     default: generated.bam
     inputBinding:
       prefix: --OUTPUT
       position: 10
   - id: addMateCigar
     label: addMateCigar
     doc: Adds the mate CIGAR tag (MC)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ADD_MATE_CIGAR
       position: 11
   - id: alignedReadsOnly
     label: alignedReadsOnly
     doc: Whether to output only aligned reads.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ALIGNED_READS_ONLY
       position: 11
   - id: alignerProperPairFlags
     label: alignerProperPairFlags
     doc: |-
       Use the aligner's idea of what a proper pair is rather than computing in this program.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ALIGNER_PROPER_PAIR_FLAGS
       position: 11
   - id: argumentsFile
     label: argumentsFile
     doc: read one or more arguments files and add them to the command line
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       prefix: --arguments_file
       position: 11
   - id: attributesToRemove
     label: attributesToRemove
     doc: Attributes from the alignment record that should be removed when merging.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --ATTRIBUTES_TO_REMOVE
       position: 11
   - id: attributesToRetain
     label: attributesToRetain
     doc: |-
       Reserved alignment attributes (tags starting with X, Y, or Z) that should be brought over from the alignment data when merging.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --ATTRIBUTES_TO_RETAIN
       position: 11
   - id: attributesToReverse
     label: attributesToReverse
     doc: Attributes on negative strand reads that need to be reversed.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --ATTRIBUTES_TO_REVERSE
       position: 11
   - id: attributesToReverseComplement
     label: attributesToReverseComplement
     doc: Attributes on negative strand reads that need to be reverse complemented.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --ATTRIBUTES_TO_REVERSE_COMPLEMENT
       position: 11
   - id: clipAdapter
     label: clipAdapter
     doc: Whether to clip adapters where identified.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --CLIP_ADAPTERS
       position: 11
   - id: clipOverlappingReads
     label: clipOverlappingReads
     doc: |-
       For paired reads, soft clip the 3' end of each read if necessary so that it does not extend past the 5' end of its mate.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --CLIP_OVERLAPPING_READS
       position: 11
   - id: expectedOrientations
     label: expectedOrientations
     doc: The expected orientation of proper read pairs.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --EXPECTED_ORIENTATIONS
       position: 11
   - id: includeSecondaryAlginments
     label: includeSecondaryAlginments
     doc: If false, do not write secondary alignments to output.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --INCLUDE_SECONDARY_ALIGNMENTS
       position: 11
   - id: isBisulfiteSequencing
     label: isBisulfiteSequencing
     doc: Whether the lane is bisulfite sequence (used when calculating the NM tag).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --IS_BISULFITE_SEQUENCE
       position: 11
   - id: matchingDictionaryTags
     label: matchingDictionaryTags
     doc: |-
       List of Sequence Records tags that must be equal (if present) in the reference dictionary and in the aligned file.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --MATCHING_DICTIONARY_TAGS
       position: 11
   - id: maxInsertionsOrDeletions
     label: maxInsertionsOrDeletions
     doc: |-
       The maximum number of insertions or deletions permitted for an alignment to be included.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --MAX_INSERTIONS_OR_DELETIONS
       position: 11
   - id: minUnclippedBases
     label: minUnclippedBases
     doc: |-
       If UNMAP_CONTAMINANT_READS is set, require this many unclipped bases or else the read will be marked as contaminant.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --MIN_UNCLIPPED_BASES
       position: 11
   - id: primaryAlignmentStrategy
     label: primaryAlignmentStrategy
     doc: |-
       Strategy for selecting primary alignment when the aligner has provided more than one alignment for a pair or fragment, and none are marked as primary, more than one is marked as primary, or the primary alignment is filtered out for some reason.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --PRIMARY_ALIGNMENT_STRATEGY
       position: 11
   - id: programGroupCommandLine
     label: programGroupCommandLine
     doc: The command line of the program group.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --PROGRAM_GROUP_COMMAND_LINE
       position: 11
   - id: programGroupName
     label: programGroupName
     doc: The name of the program group.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --PROGRAM_GROUP_NAME
       position: 11
   - id: programGroupVersion
     label: programGroupVersion
     doc: The version of the program group.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --PROGRAM_GROUP_VERSION
       position: 11
   - id: programRecordId
     label: programRecordId
     doc: The program group ID of the aligner.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --PROGRAM_RECORD_ID
       position: 11
   - id: sortOrder
     label: sortOrder
     doc: |-
       The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -SO
       position: 10
   - id: unmapContaminantReads
     label: unmapContaminantReads
     doc: |-
       Detect reads originating from foreign organisms (e.g. bacterial DNA in a non-bacterial sample),and unmap + label those reads accordingly.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --UNMAP_CONTAMINANT_READS
       position: 11
   - id: unmappedReadStrategy
     label: unmappedReadStrategy
     doc: |-
       How to deal with alignment information in reads that are being unmapped (e.g. due to cross-species contamination.) Currently ignored unless UNMAP_CONTAMINANT_READS = true.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --UNMAPPED_READ_STRATEGY
       position: 11
   - id: addPgTagToReads
     label: addPgTagToReads
     doc: Add PG tag to each read in a SAM or BAM
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ADD_PG_TAG_TO_READS
       position: 11
   - id: compressionLevel
     label: compressionLevel
     doc: Compression level for all compressed files created (e.g. BAM and GELI).
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --COMPRESSION_LEVEL
       position: 11
   - id: createIndex
     label: createIndex
     doc: Whether to create a BAM index when writing a coordinate-sorted BAM file.
     type: boolean
     default: true
     inputBinding:
       prefix: --CREATE_INDEX
       position: 11
   - id: createMd5File
     label: createMd5File
     doc: Whether to create an MD5 digest for any BAM or FASTQ files created.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --CREATE_MD5_FILE
       position: 11
   - id: maxRecordsInRam
     label: maxRecordsInRam
     doc: |-
       When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --MAX_RECORDS_IN_RAM
       position: 11
   - id: quiet
     label: quiet
     doc: Whether to suppress job-summary info on System.err.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --QUIET
       position: 11
   - id: tmpDir
     label: tmpDir
     doc: Undocumented option
     type: string
     default: /tmp/
     inputBinding:
       prefix: --TMP_DIR
       position: 11
   - id: useJdkDeflater
     label: useJdkDeflater
     doc: Whether to use the JdkDeflater (as opposed to IntelDeflater)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --use_jdk_deflater
       position: 11
   - id: useJdkInflater
     label: useJdkInflater
     doc: Whether to use the JdkInflater (as opposed to IntelInflater)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --use_jdk_inflater
       position: 11
   - id: validationStringency
     label: validationStringency
     doc: |-
       Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --VALIDATION_STRINGENCY
       position: 11
   - id: verbosity
     label: verbosity
     doc: |-
       The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --verbosity
       position: 11

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.bam
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - MergeBamAlignment
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 4, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4MergeBamAlignment


