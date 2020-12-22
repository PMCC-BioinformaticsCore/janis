:orphan:

GATK4: GatherBamFiles
===========================================

``Gatk4GatherBamFiles`` · *1 contributor · 3 versions*

b'USAGE: GatherBamFiles [arguments]
<p>Concatenate efficiently BAM files that resulted from a scattered parallel analysis.</p><p>This tool performs a rapid
'gather' or concatenation on BAM files. This is often needed in operations that have been run in parallel across
genomics regions by scattering their execution across computing nodes and cores thus resulting in smaller BAM
files.</p><p>This tool does not support SAM files</p><h3>Inputs</h3><p>A list of BAM files to combine using the INPUT
argument. These files must be provided in the order that they should be concatenated.</p><h3>Output</h3><p>A single BAM
file. The header is copied from the first input file.</p><h3>Usage example:</h3><pre>java -jar picard.jar GatherBamFiles
\
I=input1.bam \
I=input2.bam \
O=gathered_files.bam</pre><h3>Notes</h3><p>Operates via copying of the gzip blocks directly for speed but also supports
generation of an MD5 on the output and indexing of the output BAM file.</p><hr/>
Version:4.1.3.0



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.gatherbamfiles.versions import Gatk4GatherBamFiles_4_1_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4gatherbamfiles_step",
           Gatk4GatherBamFiles_4_1_3(

           )
       )
       wf.output("out", source=gatk4gatherbamfiles_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4GatherBamFiles:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4GatherBamFiles > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run Gatk4GatherBamFiles with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4GatherBamFiles





Information
------------

:ID: ``Gatk4GatherBamFiles``
:URL: *No URL to the documentation was provided*
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.3.0
:Authors: Michael Franklin
:Citations: None
:Created: 2020-05-18
:Updated: 2020-05-18


Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedBam
======  ==========  ===============


Additional configuration (inputs)
---------------------------------

=====================  =======================  =======================  ==========  ============================================================================================================================================================================================================================================================================================================
name                   type                     prefix                   position    documentation
=====================  =======================  =======================  ==========  ============================================================================================================================================================================================================================================================================================================
javaOptions            Optional<Array<String>>
compression_level      Optional<Integer>                                             Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
bams                   Optional<Array<BAM>>     --INPUT                              (-I) Two or more BAM files or text files containing lists of BAM files (one per line). This argument must be specified at least once. Required.
outputFilename         Optional<Filename>       --OUTPUT                             (-O) The output BAM file to write to. Required.
arguments_file         Optional<File>           --arguments_file                     read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null.
create_index           Optional<Boolean>        --CREATE_INDEX                       Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value: false. Possible values: {true, false}
create_md5_file        Optional<Boolean>        --CREATE_MD5_FILE                    Whether to create an MD5 digest for any BAM or FASTQ files created. Default value: false. Possible values: {true, false}
ga4gh_client_secrets   Optional<Boolean>        --GA4GH_CLIENT_SECRETS               Default value: client_secrets.json.
help                   Optional<Boolean>        --help                               (-h) display the help message Default value: false. Possible values: {true, false}
max_records_in_ram     Optional<Integer>        --MAX_RECORDS_IN_RAM                 When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.  Default value: 500000.
quiet                  Optional<Boolean>        --QUIET                              Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false}
reference_sequence     Optional<File>           --REFERENCE_SEQUENCE                 (-R) Reference sequence file. Default value: null.
tmp_dir                Optional<File>           --TMP_DIR                            One or more directories with space available to be used by this program for temporary storage of working files  This argument may be specified 0 or more times. Default value: null.
use_jdk_deflater       Optional<Boolean>        --USE_JDK_DEFLATER                   (-use_jdk_deflater)  Use the JDK Deflater instead of the Intel Deflater for writing compressed output  Default value: false. Possible values: {true, false}
use_jdk_inflater       Optional<Boolean>        --USE_JDK_INFLATER                   (-use_jdk_inflater)  Use the JDK Inflater instead of the Intel Inflater for reading compressed input  Default value: false. Possible values: {true, false}
validation_stringency  Optional<Boolean>        --VALIDATION_STRINGENCY              Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. Possible values: {STRICT, LENIENT, SILENT}
verbosity              Optional<Boolean>        --VERBOSITY                          Control verbosity of logging. Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG}
version                Optional<Boolean>        --version                            display the version number for this tool Default value: false. Possible values: {true, false}
showhidden             Optional<Boolean>        --showHidden                         (-showHidden)  display hidden arguments  Default value: false. Possible values: {true, false}
=====================  =======================  =======================  ==========  ============================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4GatherBamFiles {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       Array[File]? bams
       String? outputFilename
       File? arguments_file
       Boolean? create_index
       Boolean? create_md5_file
       Boolean? ga4gh_client_secrets
       Boolean? help
       Int? max_records_in_ram
       Boolean? quiet
       File? reference_sequence
       File? tmp_dir
       Boolean? use_jdk_deflater
       Boolean? use_jdk_inflater
       Boolean? validation_stringency
       Boolean? verbosity
       Boolean? version
       Boolean? showhidden
     }
     command <<<
       set -e
       gatk GatherBamFiles \
         --java-options '-Xmx~{((select_first([runtime_memory, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         ~{if (defined(bams) && length(select_first([bams])) > 0) then "--INPUT '" + sep("' --INPUT '", select_first([bams])) + "'" else ""} \
         --OUTPUT '~{select_first([outputFilename, "generated.merged.bam"])}' \
         ~{if defined(arguments_file) then ("--arguments_file '" + arguments_file + "'") else ""} \
         ~{if select_first([create_index, true]) then "--CREATE_INDEX" else ""} \
         ~{if (defined(create_md5_file) && select_first([create_md5_file])) then "--CREATE_MD5_FILE" else ""} \
         ~{if (defined(ga4gh_client_secrets) && select_first([ga4gh_client_secrets])) then "--GA4GH_CLIENT_SECRETS" else ""} \
         ~{if (defined(help) && select_first([help])) then "--help" else ""} \
         ~{if defined(max_records_in_ram) then ("--MAX_RECORDS_IN_RAM " + max_records_in_ram) else ''} \
         ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
         ~{if defined(reference_sequence) then ("--REFERENCE_SEQUENCE '" + reference_sequence + "'") else ""} \
         ~{if defined(tmp_dir) then ("--TMP_DIR '" + tmp_dir + "'") else ""} \
         ~{if (defined(use_jdk_deflater) && select_first([use_jdk_deflater])) then "--USE_JDK_DEFLATER" else ""} \
         ~{if (defined(use_jdk_inflater) && select_first([use_jdk_inflater])) then "--USE_JDK_INFLATER" else ""} \
         ~{if (defined(validation_stringency) && select_first([validation_stringency])) then "--VALIDATION_STRINGENCY" else ""} \
         ~{if (defined(verbosity) && select_first([verbosity])) then "--VERBOSITY" else ""} \
         ~{if (defined(version) && select_first([version])) then "--version" else ""} \
         ~{if (defined(showhidden) && select_first([showhidden])) then "--showHidden" else ""}
       if [ -f $(echo '~{select_first([outputFilename, "generated.merged.bam"])}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{select_first([outputFilename, "generated.merged.bam"])}' | sed 's/\.[^.]*$//').bai $(echo '~{select_first([outputFilename, "generated.merged.bam"])}' ).bai; fi
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.3.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.merged.bam"])
       File out_bai = select_first([outputFilename, "generated.merged.bam"]) + ".bai"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: GatherBamFiles'
   doc: |
     b'USAGE: GatherBamFiles [arguments]
     <p>Concatenate efficiently BAM files that resulted from a scattered parallel analysis.</p><p>This tool performs a rapid
     'gather' or concatenation on BAM files. This is often needed in operations that have been run in parallel across
     genomics regions by scattering their execution across computing nodes and cores thus resulting in smaller BAM
     files.</p><p>This tool does not support SAM files</p><h3>Inputs</h3><p>A list of BAM files to combine using the INPUT
     argument. These files must be provided in the order that they should be concatenated.</p><h3>Output</h3><p>A single BAM
     file. The header is copied from the first input file.</p><h3>Usage example:</h3><pre>java -jar picard.jar GatherBamFiles
     \
     I=input1.bam \
     I=input2.bam \
     O=gathered_files.bam</pre><h3>Notes</h3><p>Operates via copying of the gzip blocks directly for speed but also supports
     generation of an MD5 on the output and indexing of the output BAM file.</p><hr/>
     Version:4.1.3.0

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.3.0

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
   - id: bams
     label: bams
     doc: |-
       (-I) Two or more BAM files or text files containing lists of BAM files (one per line). This argument must be specified at least once. Required. 
     type:
     - type: array
       inputBinding:
         prefix: --INPUT
         separate: true
       items: File
     - 'null'
     inputBinding: {}
   - id: outputFilename
     label: outputFilename
     doc: (-O) The output BAM file to write to. Required.
     type:
     - string
     - 'null'
     default: generated.merged.bam
     inputBinding:
       prefix: --OUTPUT
       separate: true
   - id: arguments_file
     label: arguments_file
     doc: |-
       read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null. 
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --arguments_file
       separate: true
   - id: create_index
     label: create_index
     doc: |-
       Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value: false. Possible values: {true, false} 
     type: boolean
     default: true
     inputBinding:
       prefix: --CREATE_INDEX
       separate: true
   - id: create_md5_file
     label: create_md5_file
     doc: |-
       Whether to create an MD5 digest for any BAM or FASTQ files created. Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --CREATE_MD5_FILE
       separate: true
   - id: ga4gh_client_secrets
     label: ga4gh_client_secrets
     doc: 'Default value: client_secrets.json.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --GA4GH_CLIENT_SECRETS
       separate: true
   - id: help
     label: help
     doc: |-
       (-h) display the help message Default value: false. Possible values: {true, false}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --help
       separate: true
   - id: max_records_in_ram
     label: max_records_in_ram
     doc: |-
       When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.  Default value: 500000. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --MAX_RECORDS_IN_RAM
       separate: true
   - id: quiet
     label: quiet
     doc: |-
       Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --QUIET
       separate: true
   - id: reference_sequence
     label: reference_sequence
     doc: '(-R) Reference sequence file. Default value: null.'
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --REFERENCE_SEQUENCE
       separate: true
   - id: tmp_dir
     label: tmp_dir
     doc: |-
       One or more directories with space available to be used by this program for temporary storage of working files  This argument may be specified 0 or more times. Default value: null. 
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --TMP_DIR
       separate: true
   - id: use_jdk_deflater
     label: use_jdk_deflater
     doc: |-
       (-use_jdk_deflater)  Use the JDK Deflater instead of the Intel Deflater for writing compressed output  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --USE_JDK_DEFLATER
       separate: true
   - id: use_jdk_inflater
     label: use_jdk_inflater
     doc: |-
       (-use_jdk_inflater)  Use the JDK Inflater instead of the Intel Inflater for reading compressed input  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --USE_JDK_INFLATER
       separate: true
   - id: validation_stringency
     label: validation_stringency
     doc: |2-
        Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. Possible values: {STRICT, LENIENT, SILENT} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --VALIDATION_STRINGENCY
       separate: true
   - id: verbosity
     label: verbosity
     doc: |-
       Control verbosity of logging. Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --VERBOSITY
       separate: true
   - id: version
     label: version
     doc: |-
       display the version number for this tool Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --version
       separate: true
   - id: showhidden
     label: showhidden
     doc: |-
       (-showHidden)  display hidden arguments  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --showHidden
       separate: true

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
       glob: generated.merged.bam
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - GatherBamFiles
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4GatherBamFiles


