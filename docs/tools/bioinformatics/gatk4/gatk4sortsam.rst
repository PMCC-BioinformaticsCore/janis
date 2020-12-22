:orphan:

GATK4: SortSAM
=============================

``Gatk4SortSam`` · *1 contributor · 4 versions*

Sorts a SAM/BAM/CRAM file.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.sortsam.versions import Gatk4SortSam_4_1_4

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4sortsam_step",
           Gatk4SortSam_4_1_4(
               bam=None,
               sortOrder=None,
           )
       )
       wf.output("out", source=gatk4sortsam_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4SortSam:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4SortSam > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       sortOrder: <value>




5. Run Gatk4SortSam with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4SortSam





Information
------------

:ID: ``Gatk4SortSam``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.3/org_broadinstitute_hellbender_tools_picard_sam_SortSam.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.3/org_broadinstitute_hellbender_tools_picard_sam_SortSam.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.4.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24


Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedBam
======  ==========  ===============


Additional configuration (inputs)
---------------------------------

====================  ==========================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
name                  type                        prefix                     position  documentation
====================  ==========================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
bam                   BAM                         -I                               10  The SAM/BAM/CRAM file to sort.
sortOrder             String                      -SO                              10  The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
javaOptions           Optional<Array<String>>
compression_level     Optional<Integer>                                                Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
outputFilename        Optional<Filename>          -O                               10  The sorted SAM/BAM/CRAM output file.
argumentsFile         Optional<Array<File>>       --arguments_file                 10  read one or more arguments files and add them to the command line
compressionLevel      Optional<Integer>           --COMPRESSION_LEVEL              11  Compression level for all compressed files created (e.g. BAM and GELI).
createIndex           Optional<Boolean>           --CREATE_INDEX                   11  Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File         Optional<Boolean>           --CREATE_MD5_FILE                11  Whether to create an MD5 digest for any BAM or FASTQ files created.
maxRecordsInRam       Optional<Integer>           --MAX_RECORDS_IN_RAM             11  When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
quiet                 Optional<Boolean>           --QUIET                          11  Whether to suppress job-summary info on System.err.
reference             Optional<FastaWithIndexes>  --reference                      11  Reference sequence file.
tmpDir                Optional<String>            --TMP_DIR                        11  Undocumented option
useJdkDeflater        Optional<Boolean>           --use_jdk_deflater               11  Whether to use the JdkDeflater (as opposed to IntelDeflater)
useJdkInflater        Optional<Boolean>           --use_jdk_inflater               11  Whether to use the JdkInflater (as opposed to IntelInflater)
validationStringency  Optional<String>            --VALIDATION_STRINGENCY          11  Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
verbosity             Optional<String>            --verbosity                      11  The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
====================  ==========================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4SortSam {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       File bam
       String? outputFilename
       String sortOrder
       Array[File]? argumentsFile
       Int? compressionLevel
       Boolean? createIndex
       Boolean? createMd5File
       Int? maxRecordsInRam
       Boolean? quiet
       File? reference
       File? reference_fai
       File? reference_amb
       File? reference_ann
       File? reference_bwt
       File? reference_pac
       File? reference_sa
       File? reference_dict
       String? tmpDir
       Boolean? useJdkDeflater
       Boolean? useJdkInflater
       String? validationStringency
       String? verbosity
     }
     command <<<
       set -e
       gatk SortSam \
         --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         -I '~{bam}' \
         -O '~{select_first([outputFilename, "~{basename(bam, ".bam")}.sorted.bam"])}' \
         -SO '~{sortOrder}' \
         ~{if (defined(argumentsFile) && length(select_first([argumentsFile])) > 0) then "--arguments_file '" + sep("' --arguments_file '", select_first([argumentsFile])) + "'" else ""} \
         ~{if defined(compressionLevel) then ("--COMPRESSION_LEVEL " + compressionLevel) else ''} \
         ~{if select_first([createIndex, true]) then "--CREATE_INDEX" else ""} \
         ~{if (defined(createMd5File) && select_first([createMd5File])) then "--CREATE_MD5_FILE" else ""} \
         ~{if defined(maxRecordsInRam) then ("--MAX_RECORDS_IN_RAM " + maxRecordsInRam) else ''} \
         ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
         ~{if defined(reference) then ("--reference '" + reference + "'") else ""} \
         ~{if defined(select_first([tmpDir, "/tmp/"])) then ("--TMP_DIR '" + select_first([tmpDir, "/tmp/"]) + "'") else ""} \
         ~{if (defined(useJdkDeflater) && select_first([useJdkDeflater])) then "--use_jdk_deflater" else ""} \
         ~{if (defined(useJdkInflater) && select_first([useJdkInflater])) then "--use_jdk_inflater" else ""} \
         ~{if defined(validationStringency) then ("--VALIDATION_STRINGENCY '" + validationStringency + "'") else ""} \
         ~{if defined(verbosity) then ("--verbosity '" + verbosity + "'") else ""}
       if [ -f $(echo '~{select_first([outputFilename, "~{basename(bam, ".bam")}.sorted.bam"])}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{select_first([outputFilename, "~{basename(bam, ".bam")}.sorted.bam"])}' | sed 's/\.[^.]*$//').bai $(echo '~{select_first([outputFilename, "~{basename(bam, ".bam")}.sorted.bam"])}' ).bai; fi
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.4.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "~{basename(bam, ".bam")}.sorted.bam"])
       File out_bai = select_first([outputFilename, "~{basename(bam, ".bam")}.sorted.bam"]) + ".bai"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: SortSAM'
   doc: Sorts a SAM/BAM/CRAM file.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.4.0

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
   - id: bam
     label: bam
     doc: The SAM/BAM/CRAM file to sort.
     type: File
     inputBinding:
       prefix: -I
       position: 10
   - id: outputFilename
     label: outputFilename
     doc: The sorted SAM/BAM/CRAM output file.
     type:
     - string
     - 'null'
     default: generated.sorted.bam
     inputBinding:
       prefix: -O
       position: 10
       valueFrom: $(inputs.bam.basename.replace(/.bam$/, "")).sorted.bam
   - id: sortOrder
     label: sortOrder
     doc: |-
       The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
     type: string
     inputBinding:
       prefix: -SO
       position: 10
   - id: argumentsFile
     label: argumentsFile
     doc: read one or more arguments files and add them to the command line
     type:
     - type: array
       inputBinding:
         prefix: --arguments_file
       items: File
     - 'null'
     inputBinding:
       position: 10
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
       prefix: --reference
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
     secondaryFiles:
     - |-
       ${

               function resolveSecondary(base, secPattern) {
                 if (secPattern[0] == "^") {
                   var spl = base.split(".");
                   var endIndex = spl.length > 1 ? spl.length - 1 : 1;
                   return resolveSecondary(spl.slice(undefined, endIndex).join("."), secPattern.slice(1));
                 }
                 return base + secPattern
               }
               return [
                       {
                           path: resolveSecondary(self.path, "^.bai"),
                           basename: resolveSecondary(self.basename, ".bai"),
                           class: "File",
                       }
               ];

       }
     outputBinding:
       glob: $(inputs.bam.basename.replace(/.bam$/, "")).sorted.bam
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - SortSam
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4SortSam


