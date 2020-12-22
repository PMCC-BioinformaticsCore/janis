:orphan:

GATK4: Gather VCFs
====================================

``Gatk4GatherVcfs`` · *1 contributor · 4 versions*

GatherVcfs (Picard)
            
Gathers multiple VCF files from a scatter operation into a single VCF file. 
Input files must be supplied in genomic order and must not have events at overlapping positions.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.gathervcfs.versions import Gatk4GatherVcfs_4_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4gathervcfs_step",
           Gatk4GatherVcfs_4_1_2(
               vcfs=None,
           )
       )
       wf.output("out", source=gatk4gathervcfs_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4GatherVcfs:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4GatherVcfs > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcfs:
       - vcfs_0.vcf
       - vcfs_1.vcf




5. Run Gatk4GatherVcfs with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4GatherVcfs





Information
------------

:ID: ``Gatk4GatherVcfs``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.12.0/picard_vcf_GatherVcfs.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.12.0/picard_vcf_GatherVcfs.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.2.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-05-01
:Updated: 2019-05-01


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============


Additional configuration (inputs)
---------------------------------

====================  =======================  =======================  ==========  ======================================================================================================================================================================================================================================================================
name                  type                     prefix                   position    documentation
====================  =======================  =======================  ==========  ======================================================================================================================================================================================================================================================================
vcfs                  Array<VCF>               --INPUT                              [default: []] (-I) Input VCF file(s).
javaOptions           Optional<Array<String>>
compression_level     Optional<Integer>                                             Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
outputFilename        Optional<Filename>       --OUTPUT                             [default: null] (-O) Output VCF file.
argumentsFile         Optional<Array<File>>    --arguments_file                     [default: []] read one or more arguments files and add them to the command line
compressionLevel      Optional<Integer>        --COMPRESSION_LEVEL                  [default: 5] Compression level for all compressed files created (e.g. BAM and VCF).
createIndex           Optional<Boolean>        --CREATE_INDEX                       [default: TRUE] Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File         Optional<Boolean>        --CREATE_MD5_FILE                    [default: FALSE] Whether to create an MD5 digest for any BAM or FASTQ files created.
ga4ghClientSecrets    Optional<File>           --GA4GH_CLIENT_SECRETS               [default: client_secrets.json] Google Genomics API client_secrets.json file path.
maxRecordsInRam       Optional<Integer>        --MAX_RECORDS_IN_RAM                 [default: 500000] When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.
quiet                 Optional<Boolean>        --QUIET                              [default: FALSE] Whether to suppress job-summary info on System.err.
referenceSequence     Optional<File>           --REFERENCE_SEQUENCE                 [default: null] Reference sequence file.
tmpDir                Optional<String>         --TMP_DIR                            [default: []] One or more directories with space available to be used by this program for temporary storage of working files
useJdkDeflater        Optional<Boolean>        --USE_JDK_DEFLATER                   [default: FALSE] (-use_jdk_deflater) Use the JDK Deflater instead of the Intel Deflater for writing compressed output
useJdkInflater        Optional<Boolean>        --USE_JDK_INFLATER                   [default: FALSE] (-use_jdk_inflater) Use the JDK Inflater instead of the Intel Inflater for reading compressed input
validationStringency  Optional<String>         --VALIDATION_STRINGENCY              [default: STRICT] Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
verbosity             Optional<Boolean>        --VERBOSITY                          [default: INFO] Control verbosity of logging.
====================  =======================  =======================  ==========  ======================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4GatherVcfs {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       Array[File] vcfs
       String? outputFilename
       Array[File]? argumentsFile
       Int? compressionLevel
       Boolean? createIndex
       Boolean? createMd5File
       File? ga4ghClientSecrets
       Int? maxRecordsInRam
       Boolean? quiet
       File? referenceSequence
       String? tmpDir
       Boolean? useJdkDeflater
       Boolean? useJdkInflater
       String? validationStringency
       Boolean? verbosity
     }
     command <<<
       set -e
       gatk GatherVcfs \
         --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         ~{if length(vcfs) > 0 then "--INPUT '" + sep("' --INPUT '", vcfs) + "'" else ""} \
         --OUTPUT '~{select_first([outputFilename, "generated.gathered.vcf"])}' \
         ~{if (defined(argumentsFile) && length(select_first([argumentsFile])) > 0) then "--arguments_file '" + sep("' '", select_first([argumentsFile])) + "'" else ""} \
         ~{if defined(compressionLevel) then ("--COMPRESSION_LEVEL " + compressionLevel) else ''} \
         ~{if (defined(createIndex) && select_first([createIndex])) then "--CREATE_INDEX" else ""} \
         ~{if (defined(createMd5File) && select_first([createMd5File])) then "--CREATE_MD5_FILE" else ""} \
         ~{if defined(ga4ghClientSecrets) then ("--GA4GH_CLIENT_SECRETS '" + ga4ghClientSecrets + "'") else ""} \
         ~{if defined(maxRecordsInRam) then ("--MAX_RECORDS_IN_RAM " + maxRecordsInRam) else ''} \
         ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
         ~{if defined(referenceSequence) then ("--REFERENCE_SEQUENCE '" + referenceSequence + "'") else ""} \
         ~{if defined(select_first([tmpDir, "/tmp"])) then ("--TMP_DIR '" + select_first([tmpDir, "/tmp"]) + "'") else ""} \
         ~{if (defined(useJdkDeflater) && select_first([useJdkDeflater])) then "--USE_JDK_DEFLATER" else ""} \
         ~{if (defined(useJdkInflater) && select_first([useJdkInflater])) then "--USE_JDK_INFLATER" else ""} \
         ~{if defined(validationStringency) then ("--VALIDATION_STRINGENCY '" + validationStringency + "'") else ""} \
         ~{if (defined(verbosity) && select_first([verbosity])) then "--VERBOSITY" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.2.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.gathered.vcf"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: Gather VCFs'
   doc: |-
     GatherVcfs (Picard)
              
     Gathers multiple VCF files from a scatter operation into a single VCF file. 
     Input files must be supplied in genomic order and must not have events at overlapping positions.

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
   - id: vcfs
     label: vcfs
     doc: '[default: []] (-I) Input VCF file(s).'
     type:
       type: array
       inputBinding:
         prefix: --INPUT
       items: File
     inputBinding: {}
   - id: outputFilename
     label: outputFilename
     doc: '[default: null] (-O) Output VCF file.'
     type:
     - string
     - 'null'
     default: generated.gathered.vcf
     inputBinding:
       prefix: --OUTPUT
   - id: argumentsFile
     label: argumentsFile
     doc: '[default: []] read one or more arguments files and add them to the command
       line'
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       prefix: --arguments_file
   - id: compressionLevel
     label: compressionLevel
     doc: |-
       [default: 5] Compression level for all compressed files created (e.g. BAM and VCF).
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --COMPRESSION_LEVEL
   - id: createIndex
     label: createIndex
     doc: |-
       [default: TRUE] Whether to create a BAM index when writing a coordinate-sorted BAM file.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --CREATE_INDEX
   - id: createMd5File
     label: createMd5File
     doc: |-
       [default: FALSE] Whether to create an MD5 digest for any BAM or FASTQ files created.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --CREATE_MD5_FILE
   - id: ga4ghClientSecrets
     label: ga4ghClientSecrets
     doc: |-
       [default: client_secrets.json] Google Genomics API client_secrets.json file path.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --GA4GH_CLIENT_SECRETS
   - id: maxRecordsInRam
     label: maxRecordsInRam
     doc: |-
       [default: 500000] When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --MAX_RECORDS_IN_RAM
   - id: quiet
     label: quiet
     doc: '[default: FALSE] Whether to suppress job-summary info on System.err.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --QUIET
   - id: referenceSequence
     label: referenceSequence
     doc: '[default: null] Reference sequence file.'
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --REFERENCE_SEQUENCE
   - id: tmpDir
     label: tmpDir
     doc: |-
       [default: []] One or more directories with space available to be used by this program for temporary storage of working files
     type: string
     default: /tmp
     inputBinding:
       prefix: --TMP_DIR
   - id: useJdkDeflater
     label: useJdkDeflater
     doc: |-
       [default: FALSE] (-use_jdk_deflater) Use the JDK Deflater instead of the Intel Deflater for writing compressed output
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --USE_JDK_DEFLATER
   - id: useJdkInflater
     label: useJdkInflater
     doc: |-
       [default: FALSE] (-use_jdk_inflater) Use the JDK Inflater instead of the Intel Inflater for reading compressed input
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --USE_JDK_INFLATER
   - id: validationStringency
     label: validationStringency
     doc: |-
       [default: STRICT] Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --VALIDATION_STRINGENCY
   - id: verbosity
     label: verbosity
     doc: '[default: INFO] Control verbosity of logging.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --VERBOSITY

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.gathered.vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - GatherVcfs
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4GatherVcfs


