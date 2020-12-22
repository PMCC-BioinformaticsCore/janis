:orphan:

GATK4: Convert a FASTQ file to an unaligned BAM or SAM file.
==============================================================================

``Gatk4FastqToSam`` · *2 contributors · 4 versions*

Converts a FASTQ file to an unaligned BAM or SAM file.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.fastqtosam.versions import Gatk4FastqToSam_4_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4fastqtosam_step",
           Gatk4FastqToSam_4_1_2(
               fastqR1=None,
           )
       )
       wf.output("out", source=gatk4fastqtosam_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4FastqToSam:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4FastqToSam > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fastqR1: fastqR1.fastq.gz




5. Run Gatk4FastqToSam with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4FastqToSam





Information
------------

:ID: ``Gatk4FastqToSam``
:URL: `https://gatk.broadinstitute.org/hc/en-us/articles/360037226792-FastqToSam-Picard- <https://gatk.broadinstitute.org/hc/en-us/articles/360037226792-FastqToSam-Picard->`_
:Versions: 4.1.4.1, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.2.0
:Authors: Michael Franklin (@illusional), Matthias De Smet(@matthdsm)
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2020-02-26
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

========================  ==========================  ==============================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
name                      type                        prefix                            position  documentation
========================  ==========================  ==============================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
fastqR1                   FastqGz                     --FASTQ                                 10  Input fastq file (optionally gzipped) for single end data, or first read in paired end data.
javaOptions               Optional<Array<String>>
compression_level         Optional<Integer>                                                       Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
fastqR2                   Optional<FastqGz>           --FASTQ2                                10  Input fastq file (optionally gzipped) for single end data, or first read in paired end data.
sampleName                Optional<String>            --SAMPLE_NAME                           10  Input fastq file (optionally gzipped) for single end data, or first read in paired end data.
reference                 Optional<FastaWithIndexes>  --REFERENCE_SEQUENCE                    10  Reference sequence file.
outputFilename            Optional<Filename>          --OUTPUT                                10  Merged SAM or BAM file to write to.
allowAndIgnoreEmptyLines  Optional<Boolean>           --ALLOW_AND_IGNORE_EMPTY_LINES          11  Allow (and ignore) empty lines
argumentsFile             Optional<Array<File>>       --arguments_file                        11  read one or more arguments files and add them to the command line
comment                   Optional<Array<String>>     --COMMENT                               11  Comment(s) to include in the merged output file's header.
description               Optional<Array<String>>     --DESCRIPTION                           11  Inserted into the read group header
libraryName               Optional<Array<String>>     --LIBRARY_NAME                          11  The library name to place into the LB attribute in the read group header
maxQ                      Optional<Integer>           --MAX_Q                                 11  Maximum quality allowed in the input fastq. An exception will be thrown if a quality is greater than this value.
minQ                      Optional<Integer>           --MIN_Q                                 11  Minimum quality allowed in the input fastq. An exception will be thrown if a quality is less than this value.
platform                  Optional<String>            --PLATFORM                              11  The platform type (e.g. ILLUMINA, SOLID) to insert into the read group header.
platformModel             Optional<String>            --PLATFORM_MODEL                        11  Platform model to insert into the group header (free-form text providing further details of the platform/technology used).
platformUnit              Optional<String>            --PLATFORM_UNIT                         11  The expected orientation of proper read pairs.
predictedInsertSize       Optional<Integer>           --PREDICTED_INSERT_SIZE                 11  Predicted median insert size, to insert into the read group header.
programGroup              Optional<String>            --PROGRAM_GROUP                         11  Program group to insert into the read group header.
readGroupName             Optional<String>            --READ_GROUP_NAME                       11  Read group name.
runDate                   Optional<String>            --RUN_DATE                              11  Date the run was produced, to insert into the read group header
sequencingCenter          Optional<String>            --SEQUENCING_CENTER                     11  The sequencing center from which the data originated.
sortOrder                 Optional<String>            -SO                                     10  The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
useSequenctialFastqs      Optional<Boolean>           --USE_SEQUENTIAL_FASTQS                 11  Use sequential fastq files with the suffix _###.fastq or _###.fastq.gz.
compressionLevel          Optional<Integer>           --COMPRESSION_LEVEL                     11  Compression level for all compressed files created (e.g. BAM and GELI).
createIndex               Optional<Boolean>           --CREATE_INDEX                          11  Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File             Optional<Boolean>           --CREATE_MD5_FILE                       11  Whether to create an MD5 digest for any BAM or FASTQ files created.
maxRecordsInRam           Optional<Integer>           --MAX_RECORDS_IN_RAM                    11  When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
quiet                     Optional<Boolean>           --QUIET                                 11  Whether to suppress job-summary info on System.err.
tmpDir                    Optional<String>            --TMP_DIR                               11  Undocumented option
useJdkDeflater            Optional<Boolean>           --use_jdk_deflater                      11  Whether to use the JdkDeflater (as opposed to IntelDeflater)
useJdkInflater            Optional<Boolean>           --use_jdk_inflater                      11  Whether to use the JdkInflater (as opposed to IntelInflater)
validationStringency      Optional<String>            --VALIDATION_STRINGENCY                 11  Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
verbosity                 Optional<String>            --verbosity                             11  The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
========================  ==========================  ==============================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4FastqToSam {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       File fastqR1
       File? fastqR2
       String? sampleName
       File? reference
       File? reference_fai
       File? reference_amb
       File? reference_ann
       File? reference_bwt
       File? reference_pac
       File? reference_sa
       File? reference_dict
       String? outputFilename
       Boolean? allowAndIgnoreEmptyLines
       Array[File]? argumentsFile
       Array[String]? comment
       Array[String]? description
       Array[String]? libraryName
       Int? maxQ
       Int? minQ
       String? platform
       String? platformModel
       String? platformUnit
       Int? predictedInsertSize
       String? programGroup
       String? readGroupName
       String? runDate
       String? sequencingCenter
       String? sortOrder
       Boolean? useSequenctialFastqs
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
       gatk FastqToSam \
         --java-options '-Xmx~{((select_first([runtime_memory, 4, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         --FASTQ '~{fastqR1}' \
         ~{if defined(fastqR2) then ("--FASTQ2 '" + fastqR2 + "'") else ""} \
         ~{if defined(sampleName) then ("--SAMPLE_NAME '" + sampleName + "'") else ""} \
         ~{if defined(reference) then ("--REFERENCE_SEQUENCE '" + reference + "'") else ""} \
         --OUTPUT '~{select_first([outputFilename, "generated.bam"])}' \
         ~{if defined(sortOrder) then ("-SO '" + sortOrder + "'") else ""} \
         ~{if (defined(allowAndIgnoreEmptyLines) && select_first([allowAndIgnoreEmptyLines])) then "--ALLOW_AND_IGNORE_EMPTY_LINES" else ""} \
         ~{if (defined(argumentsFile) && length(select_first([argumentsFile])) > 0) then "--arguments_file '" + sep("' '", select_first([argumentsFile])) + "'" else ""} \
         ~{if (defined(comment) && length(select_first([comment])) > 0) then "--COMMENT '" + sep("' '", select_first([comment])) + "'" else ""} \
         ~{if (defined(description) && length(select_first([description])) > 0) then "--DESCRIPTION '" + sep("' '", select_first([description])) + "'" else ""} \
         ~{if (defined(libraryName) && length(select_first([libraryName])) > 0) then "--LIBRARY_NAME '" + sep("' '", select_first([libraryName])) + "'" else ""} \
         ~{if defined(maxQ) then ("--MAX_Q " + maxQ) else ''} \
         ~{if defined(minQ) then ("--MIN_Q " + minQ) else ''} \
         ~{if defined(platform) then ("--PLATFORM '" + platform + "'") else ""} \
         ~{if defined(platformModel) then ("--PLATFORM_MODEL '" + platformModel + "'") else ""} \
         ~{if defined(platformUnit) then ("--PLATFORM_UNIT '" + platformUnit + "'") else ""} \
         ~{if defined(predictedInsertSize) then ("--PREDICTED_INSERT_SIZE " + predictedInsertSize) else ''} \
         ~{if defined(programGroup) then ("--PROGRAM_GROUP '" + programGroup + "'") else ""} \
         ~{if defined(readGroupName) then ("--READ_GROUP_NAME '" + readGroupName + "'") else ""} \
         ~{if defined(runDate) then ("--RUN_DATE '" + runDate + "'") else ""} \
         ~{if defined(sequencingCenter) then ("--SEQUENCING_CENTER '" + sequencingCenter + "'") else ""} \
         ~{if (defined(useSequenctialFastqs) && select_first([useSequenctialFastqs])) then "--USE_SEQUENTIAL_FASTQS" else ""} \
         ~{if defined(compressionLevel) then ("--COMPRESSION_LEVEL " + compressionLevel) else ''} \
         ~{if (defined(createIndex) && select_first([createIndex])) then "--CREATE_INDEX" else ""} \
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
       docker: "broadinstitute/gatk:4.1.2.0"
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
   label: 'GATK4: Convert a FASTQ file to an unaligned BAM or SAM file.'
   doc: Converts a FASTQ file to an unaligned BAM or SAM file.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.2.0

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
   - id: fastqR1
     label: fastqR1
     doc: |-
       Input fastq file (optionally gzipped) for single end data, or first read in paired end data.
     type: File
     inputBinding:
       prefix: --FASTQ
       position: 10
   - id: fastqR2
     label: fastqR2
     doc: |-
       Input fastq file (optionally gzipped) for single end data, or first read in paired end data.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --FASTQ2
       position: 10
   - id: sampleName
     label: sampleName
     doc: |-
       Input fastq file (optionally gzipped) for single end data, or first read in paired end data.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --SAMPLE_NAME
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
   - id: allowAndIgnoreEmptyLines
     label: allowAndIgnoreEmptyLines
     doc: Allow (and ignore) empty lines
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ALLOW_AND_IGNORE_EMPTY_LINES
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
   - id: comment
     label: comment
     doc: Comment(s) to include in the merged output file's header.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --COMMENT
       position: 11
   - id: description
     label: description
     doc: Inserted into the read group header
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --DESCRIPTION
       position: 11
   - id: libraryName
     label: libraryName
     doc: The library name to place into the LB attribute in the read group header
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --LIBRARY_NAME
       position: 11
   - id: maxQ
     label: maxQ
     doc: |-
       Maximum quality allowed in the input fastq. An exception will be thrown if a quality is greater than this value.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --MAX_Q
       position: 11
   - id: minQ
     label: minQ
     doc: |-
       Minimum quality allowed in the input fastq. An exception will be thrown if a quality is less than this value.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --MIN_Q
       position: 11
   - id: platform
     label: platform
     doc: The platform type (e.g. ILLUMINA, SOLID) to insert into the read group header.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --PLATFORM
       position: 11
   - id: platformModel
     label: platformModel
     doc: |-
       Platform model to insert into the group header (free-form text providing further details of the platform/technology used).
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --PLATFORM_MODEL
       position: 11
   - id: platformUnit
     label: platformUnit
     doc: The expected orientation of proper read pairs.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --PLATFORM_UNIT
       position: 11
   - id: predictedInsertSize
     label: predictedInsertSize
     doc: Predicted median insert size, to insert into the read group header.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --PREDICTED_INSERT_SIZE
       position: 11
   - id: programGroup
     label: programGroup
     doc: Program group to insert into the read group header.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --PROGRAM_GROUP
       position: 11
   - id: readGroupName
     label: readGroupName
     doc: Read group name.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --READ_GROUP_NAME
       position: 11
   - id: runDate
     label: runDate
     doc: Date the run was produced, to insert into the read group header
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --RUN_DATE
       position: 11
   - id: sequencingCenter
     label: sequencingCenter
     doc: The sequencing center from which the data originated.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --SEQUENCING_CENTER
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
   - id: useSequenctialFastqs
     label: useSequenctialFastqs
     doc: Use sequential fastq files with the suffix _###.fastq or _###.fastq.gz.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --USE_SEQUENTIAL_FASTQS
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
     type:
     - boolean
     - 'null'
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
   - FastqToSam
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 4, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4FastqToSam


