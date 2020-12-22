:orphan:

GATK4: Mark Duplicates
============================================

``Gatk4MarkDuplicates`` · *1 contributor · 4 versions*

MarkDuplicates (Picard): Identifies duplicate reads.

This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are 
defined as originating from a single fragment of DNA. Duplicates can arise during sample 
preparation e.g. library construction using PCR. See also EstimateLibraryComplexity for 
additional notes on PCR duplication artifacts. Duplicate reads can also result from a single 
amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the 
sequencing instrument. These duplication artifacts are referred to as optical duplicates.

The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads 
and read-pairs in a SAM/BAM file. An BARCODE_TAG option is available to facilitate duplicate
marking using molecular barcodes. After duplicate reads are collected, the tool differentiates 
the primary and duplicate reads using an algorithm that ranks reads by the sums of their 
base-quality scores (default method).

The tool's main output is a new SAM or BAM file, in which duplicates have been identified 
in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, 
which corresponds to a decimal value of 1024. If you are not familiar with this type of annotation, 
please see the following blog post for additional information.

Although the bitwise flag annotation indicates whether a read was marked as a duplicate, 
it does not identify the type of duplicate. To do this, a new tag called the duplicate type (DT) 
tag was recently added as an optional output in the 'optional field' section of a SAM/BAM file. 
Invoking the TAGGING_POLICY option, you can instruct the program to mark all the duplicates (All), 
only the optical duplicates (OpticalOnly), or no duplicates (DontTag). The records within the 
output of a SAM/BAM file will have values for the 'DT' tag (depending on the invoked TAGGING_POLICY), 
as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ). 
This tool uses the READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the 
primary methods to identify and differentiate duplicate types. Set READ_NAME_REGEX to null to 
skip optical duplicate detection, e.g. for RNA-seq or other data where duplicate sets are 
extremely large and estimating library complexity is not an aim. Note that without optical 
duplicate counts, library size estimation will be inaccurate.

MarkDuplicates also produces a metrics file indicating the numbers 
of duplicates for both single- and paired-end reads.

The program can take either coordinate-sorted or query-sorted inputs, however the behavior 
is slightly different. When the input is coordinate-sorted, unmapped mates of mapped records 
and supplementary/secondary alignments are not marked as duplicates. However, when the input 
is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary 
reads are not excluded from the duplication test and can be marked as duplicate reads.

If desired, duplicates can be removed using the REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES options.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.markduplicates.versions import Gatk4MarkDuplicates_4_1_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4markduplicates_step",
           Gatk4MarkDuplicates_4_1_3(
               bam=None,
           )
       )
       wf.output("out", source=gatk4markduplicates_step.out)
       wf.output("metrics", source=gatk4markduplicates_step.metrics)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4MarkDuplicates:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4MarkDuplicates > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam:
       - bam_0.bam
       - bam_1.bam




5. Run Gatk4MarkDuplicates with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4MarkDuplicates





Information
------------

:ID: ``Gatk4MarkDuplicates``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_markduplicates_MarkDuplicates.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_markduplicates_MarkDuplicates.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.3.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24


Outputs
-----------

=======  ==========  ===============
name     type        documentation
=======  ==========  ===============
out      IndexedBam
metrics  tsv
=======  ==========  ===============


Additional configuration (inputs)
---------------------------------

=============================  =======================  ==================================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
name                           type                     prefix                                position  documentation
=============================  =======================  ==================================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
bam                            Array<BAM>               -I                                          10  One or more input SAM or BAM files to analyze. Must be coordinate sorted.
outputPrefix                   Optional<String>
outputFilename                 Optional<Filename>       -O                                          10  File to write duplication metrics to
metricsFilename                Optional<Filename>       -M                                          10  The output file to write marked records to.
javaOptions                    Optional<Array<String>>
compression_level              Optional<Integer>                                                        Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
argumentsFile                  Optional<Array<File>>    --arguments_file                            10  read one or more arguments files and add them to the command line
assumeSortOrder                Optional<String>         -ASO                                            If not null, assume that the input file has this order even if the header says otherwise. Exclusion: This argument cannot be used at the same time as ASSUME_SORTED. The --ASSUME_SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
barcodeTag                     Optional<String>         --BARCODE_TAG                                   Barcode SAM tag (ex. BC for 10X Genomics)
comment                        Optional<Array<String>>  -CO                                             Comment(s) to include in the output file's header.
createIndex                    Optional<Boolean>        --CREATE_INDEX                              11  Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File                  Optional<Boolean>        --CREATE_MD5_FILE                           11  Whether to create an MD5 digest for any BAM or FASTQ files created.
maxRecordsInRam                Optional<Integer>        --MAX_RECORDS_IN_RAM                        11  When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
quiet                          Optional<Boolean>        --QUIET                                     11  Whether to suppress job-summary info on System.err.
tmpDir                         Optional<String>         --TMP_DIR                                   11  Undocumented option
useJdkDeflater                 Optional<Boolean>        --use_jdk_deflater                          11  Whether to use the JdkDeflater (as opposed to IntelDeflater)
useJdkInflater                 Optional<Boolean>        --use_jdk_inflater                          11  Whether to use the JdkInflater (as opposed to IntelInflater)
validationStringency           Optional<String>         --VALIDATION_STRINGENCY                     11  Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
verbosity                      Optional<String>         --verbosity                                 11  The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
opticalDuplicatePixelDistance  Optional<Integer>        --OPTICAL_DUPLICATE_PIXEL_DISTANCE              The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is more appropriate. For other platforms and models, users should experiment to find what works best.
=============================  =======================  ==================================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4MarkDuplicates {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[File] bam
       String? outputPrefix
       String? outputFilename
       String? metricsFilename
       Array[String]? javaOptions
       Int? compression_level
       Array[File]? argumentsFile
       String? assumeSortOrder
       String? barcodeTag
       Array[String]? comment
       Boolean? createIndex
       Boolean? createMd5File
       Int? maxRecordsInRam
       Boolean? quiet
       String? tmpDir
       Boolean? useJdkDeflater
       Boolean? useJdkInflater
       String? validationStringency
       String? verbosity
       Int? opticalDuplicatePixelDistance
     }
     command <<<
       set -e
       gatk MarkDuplicates \
         --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         ~{if defined(assumeSortOrder) then ("-ASO '" + assumeSortOrder + "'") else ""} \
         ~{if defined(barcodeTag) then ("--BARCODE_TAG '" + barcodeTag + "'") else ""} \
         ~{if (defined(comment) && length(select_first([comment])) > 0) then "-CO '" + sep("' '", select_first([comment])) + "'" else ""} \
         ~{if defined(opticalDuplicatePixelDistance) then ("--OPTICAL_DUPLICATE_PIXEL_DISTANCE " + opticalDuplicatePixelDistance) else ''} \
         ~{if length(bam) > 0 then "-I '" + sep("' '", bam) + "'" else ""} \
         -O '~{select_first([outputFilename, "~{select_first([outputPrefix, "generated"])}.markduped.bam"])}' \
         -M '~{select_first([metricsFilename, "~{select_first([outputPrefix, "generated"])}.metrics.txt"])}' \
         ~{if (defined(argumentsFile) && length(select_first([argumentsFile])) > 0) then "--arguments_file '" + sep("' '", select_first([argumentsFile])) + "'" else ""} \
         ~{if select_first([createIndex, true]) then "--CREATE_INDEX" else ""} \
         ~{if (defined(createMd5File) && select_first([createMd5File])) then "--CREATE_MD5_FILE" else ""} \
         ~{if defined(maxRecordsInRam) then ("--MAX_RECORDS_IN_RAM " + maxRecordsInRam) else ''} \
         ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
         ~{if defined(select_first([tmpDir, "tmp/"])) then ("--TMP_DIR '" + select_first([tmpDir, "tmp/"]) + "'") else ""} \
         ~{if (defined(useJdkDeflater) && select_first([useJdkDeflater])) then "--use_jdk_deflater" else ""} \
         ~{if (defined(useJdkInflater) && select_first([useJdkInflater])) then "--use_jdk_inflater" else ""} \
         ~{if defined(validationStringency) then ("--VALIDATION_STRINGENCY '" + validationStringency + "'") else ""} \
         ~{if defined(verbosity) then ("--verbosity '" + verbosity + "'") else ""}
       if [ -f $(echo '~{select_first([outputFilename, "~{select_first([outputPrefix, "generated"])}.markduped.bam"])}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{select_first([outputFilename, "~{select_first([outputPrefix, "generated"])}.markduped.bam"])}' | sed 's/\.[^.]*$//').bai $(echo '~{select_first([outputFilename, "~{select_first([outputPrefix, "generated"])}.markduped.bam"])}' ).bai; fi
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.3.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "~{select_first([outputPrefix, "generated"])}.markduped.bam"])
       File out_bai = select_first([outputFilename, "~{select_first([outputPrefix, "generated"])}.markduped.bam"]) + ".bai"
       File metrics = select_first([metricsFilename, "~{select_first([outputPrefix, "generated"])}.metrics.txt"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: Mark Duplicates'
   doc: |-
     MarkDuplicates (Picard): Identifies duplicate reads.

     This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are 
     defined as originating from a single fragment of DNA. Duplicates can arise during sample 
     preparation e.g. library construction using PCR. See also EstimateLibraryComplexity for 
     additional notes on PCR duplication artifacts. Duplicate reads can also result from a single 
     amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the 
     sequencing instrument. These duplication artifacts are referred to as optical duplicates.

     The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads 
     and read-pairs in a SAM/BAM file. An BARCODE_TAG option is available to facilitate duplicate
     marking using molecular barcodes. After duplicate reads are collected, the tool differentiates 
     the primary and duplicate reads using an algorithm that ranks reads by the sums of their 
     base-quality scores (default method).

     The tool's main output is a new SAM or BAM file, in which duplicates have been identified 
     in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, 
     which corresponds to a decimal value of 1024. If you are not familiar with this type of annotation, 
     please see the following blog post for additional information.

     Although the bitwise flag annotation indicates whether a read was marked as a duplicate, 
     it does not identify the type of duplicate. To do this, a new tag called the duplicate type (DT) 
     tag was recently added as an optional output in the 'optional field' section of a SAM/BAM file. 
     Invoking the TAGGING_POLICY option, you can instruct the program to mark all the duplicates (All), 
     only the optical duplicates (OpticalOnly), or no duplicates (DontTag). The records within the 
     output of a SAM/BAM file will have values for the 'DT' tag (depending on the invoked TAGGING_POLICY), 
     as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ). 
     This tool uses the READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the 
     primary methods to identify and differentiate duplicate types. Set READ_NAME_REGEX to null to 
     skip optical duplicate detection, e.g. for RNA-seq or other data where duplicate sets are 
     extremely large and estimating library complexity is not an aim. Note that without optical 
     duplicate counts, library size estimation will be inaccurate.

     MarkDuplicates also produces a metrics file indicating the numbers 
     of duplicates for both single- and paired-end reads.

     The program can take either coordinate-sorted or query-sorted inputs, however the behavior 
     is slightly different. When the input is coordinate-sorted, unmapped mates of mapped records 
     and supplementary/secondary alignments are not marked as duplicates. However, when the input 
     is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary 
     reads are not excluded from the duplication test and can be marked as duplicate reads.

     If desired, duplicates can be removed using the REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES options.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.3.0

   inputs:
   - id: bam
     label: bam
     doc: One or more input SAM or BAM files to analyze. Must be coordinate sorted.
     type:
       type: array
       items: File
     inputBinding:
       prefix: -I
       position: 10
   - id: outputPrefix
     label: outputPrefix
     type:
     - string
     - 'null'
   - id: outputFilename
     label: outputFilename
     doc: File to write duplication metrics to
     type:
     - string
     - 'null'
     default: generated.markduped.bam
     inputBinding:
       prefix: -O
       position: 10
       valueFrom: |-
         $([inputs.outputPrefix, "generated"].filter(function (inner) { return inner != null })[0]).markduped.bam
   - id: metricsFilename
     label: metricsFilename
     doc: The output file to write marked records to.
     type:
     - string
     - 'null'
     default: generated.metrics.txt
     inputBinding:
       prefix: -M
       position: 10
       valueFrom: |-
         $([inputs.outputPrefix, "generated"].filter(function (inner) { return inner != null })[0]).metrics.txt
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
   - id: argumentsFile
     label: argumentsFile
     doc: read one or more arguments files and add them to the command line
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       prefix: --arguments_file
       position: 10
   - id: assumeSortOrder
     label: assumeSortOrder
     doc: |-
       If not null, assume that the input file has this order even if the header says otherwise. Exclusion: This argument cannot be used at the same time as ASSUME_SORTED. The --ASSUME_SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -ASO
   - id: barcodeTag
     label: barcodeTag
     doc: Barcode SAM tag (ex. BC for 10X Genomics)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --BARCODE_TAG
   - id: comment
     label: comment
     doc: Comment(s) to include in the output file's header.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: -CO
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
     default: tmp/
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
   - id: opticalDuplicatePixelDistance
     label: opticalDuplicatePixelDistance
     doc: |-
       The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is more appropriate. For other platforms and models, users should experiment to find what works best.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --OPTICAL_DUPLICATE_PIXEL_DISTANCE

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
       glob: |-
         $([inputs.outputPrefix, "generated"].filter(function (inner) { return inner != null })[0]).markduped.bam
       loadContents: false
   - id: metrics
     label: metrics
     type: File
     outputBinding:
       glob: |-
         $([inputs.outputPrefix, "generated"].filter(function (inner) { return inner != null })[0]).metrics.txt
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - MarkDuplicates
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4MarkDuplicates


