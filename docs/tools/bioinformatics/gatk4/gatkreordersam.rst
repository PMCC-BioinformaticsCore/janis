:orphan:

Gatk4: ReorderSam
==================================

``GatkReorderSam`` · *1 contributor · 1 version*


USAGE: ReorderSam [arguments]

Not to be confused with SortSam which sorts a SAM or BAM file with a valid sequence dictionary, 
ReorderSam reorders
reads in a SAM/BAM file to match the contig ordering in a provided reference file, 
as determined by exact name matching
of contigs.  Reads mapped to contigs absent in the new reference 
are dropped. Runs substantially faster if the input is
an indexed BAM file.

Example:

.. code-tool: none

   java -jar picard.jar ReorderSam \
       INPUT=sample.bam \
       OUTPUT=reordered.bam \
       REFERENCE=reference_with_different_order.fasta
       
Version:4.1.3.0


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.reordersam.versions import Gatk4ReorderSam_4_1_4

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatkreordersam_step",
           Gatk4ReorderSam_4_1_4(
               inp=None,
               sequence_dictionary=None,
           )
       )
       wf.output("out", source=gatkreordersam_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for GatkReorderSam:

.. code-block:: bash

   # user inputs
   janis inputs GatkReorderSam > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inp: inp.bam
       sequence_dictionary: sequence_dictionary




5. Run GatkReorderSam with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       GatkReorderSam





Information
------------

:ID: ``GatkReorderSam``
:URL: *No URL to the documentation was provided*
:Versions: 4.1.4.0
:Container: broadinstitute/gatk:4.1.4.0
:Authors: illusional
:Citations: None
:Created: 2020-05-15
:Updated: 2020-05-15


Outputs
-----------

======  ==========  ===============================
name    type        documentation
======  ==========  ===============================
out     IndexedBam  BAM to write extracted reads to
======  ==========  ===============================


Additional configuration (inputs)
---------------------------------

=================================  =======================  ===================================  ==========  ============================================================================================================================================================================================================================================================================================================
name                               type                     prefix                               position    documentation
=================================  =======================  ===================================  ==========  ============================================================================================================================================================================================================================================================================================================
inp                                BAM                      --INPUT                                          (-I) Input file (SAM or BAM) to extract reads from. Required.
sequence_dictionary                File                     --SEQUENCE_DICTIONARY                            A Sequence Dictionary for the OUTPUT file (can be read from one of the following file types (SAM, BAM, VCF, BCF, Interval List, Fasta, or Dict)
javaOptions                        Optional<Array<String>>
compression_level                  Optional<Integer>                                                         Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
outputFilename                     Optional<Filename>       --OUTPUT                                         (-O) Output file (SAM or BAM) to write extracted reads to. Required.
allow_contig_length_discordance    Optional<Boolean>        --ALLOW_CONTIG_LENGTH_DISCORDANCE                (-U)  If true, then permits mapping from a read contig to a new reference contig with the same name but a different length.  Highly dangerous, only use if you know what you are doing.  Default value: false. Possible values: {true, false}
allow_incomplete_dict_concordance  Optional<Boolean>        --ALLOW_INCOMPLETE_DICT_CONCORDANCE              (-S)  If true, then allows only a partial overlap of the original contigs with the new reference sequence contigs.  By default, this tool requires a corresponding contig in the new reference for each read contig  Default value: false. Possible values: {true, false}
arguments_file                     Optional<File>           --arguments_file                                 read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null.
create_index                       Optional<Boolean>        --CREATE_INDEX                                   Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value: false. Possible values: {true, false}
create_md5_file                    Optional<Boolean>        --CREATE_MD5_FILE                                Whether to create an MD5 digest for any BAM or FASTQ files created. Default value: false. Possible values: {true, false}
ga4gh_client_secrets               Optional<Boolean>        --GA4GH_CLIENT_SECRETS                           Default value: client_secrets.json.
help                               Optional<Boolean>        --help                                           (-h) display the help message Default value: false. Possible values: {true, false}
max_records_in_ram                 Optional<Integer>        --MAX_RECORDS_IN_RAM                             When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.  Default value: 500000.
quiet                              Optional<Boolean>        --QUIET                                          Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false}
tmp_dir                            Optional<File>           --TMP_DIR                                        One or more directories with space available to be used by this program for temporary storage of working files  This argument may be specified 0 or more times. Default value: null.
use_jdk_deflater                   Optional<Boolean>        --USE_JDK_DEFLATER                               (-use_jdk_deflater)  Use the JDK Deflater instead of the Intel Deflater for writing compressed output  Default value: false. Possible values: {true, false}
use_jdk_inflater                   Optional<Boolean>        --USE_JDK_INFLATER                               (-use_jdk_inflater)  Use the JDK Inflater instead of the Intel Inflater for reading compressed input  Default value: false. Possible values: {true, false}
validation_stringency              Optional<Boolean>        --VALIDATION_STRINGENCY                          Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. Possible values: {STRICT, LENIENT, SILENT}
verbosity                          Optional<Boolean>        --VERBOSITY                                      Control verbosity of logging. Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG}
version                            Optional<Boolean>        --version                                        display the version number for this tool Default value: false. Possible values: {true, false}
showhidden                         Optional<Boolean>        --showHidden                                     (-showHidden)  display hidden arguments  Default value: false. Possible values: {true, false}
=================================  =======================  ===================================  ==========  ============================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task GatkReorderSam {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       File inp
       String? outputFilename
       File sequence_dictionary
       Boolean? allow_contig_length_discordance
       Boolean? allow_incomplete_dict_concordance
       File? arguments_file
       Boolean? create_index
       Boolean? create_md5_file
       Boolean? ga4gh_client_secrets
       Boolean? help
       Int? max_records_in_ram
       Boolean? quiet
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
       gatk ReorderSam \
         --java-options '-Xmx~{((select_first([runtime_memory, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         --INPUT '~{inp}' \
         --OUTPUT '~{select_first([outputFilename, "~{basename(inp, ".bam")}.bam"])}' \
         --SEQUENCE_DICTIONARY '~{sequence_dictionary}' \
         ~{if (defined(allow_contig_length_discordance) && select_first([allow_contig_length_discordance])) then "--ALLOW_CONTIG_LENGTH_DISCORDANCE" else ""} \
         ~{if (defined(allow_incomplete_dict_concordance) && select_first([allow_incomplete_dict_concordance])) then "--ALLOW_INCOMPLETE_DICT_CONCORDANCE" else ""} \
         ~{if defined(arguments_file) then ("--arguments_file '" + arguments_file + "'") else ""} \
         ~{if select_first([create_index, true]) then "--CREATE_INDEX" else ""} \
         ~{if (defined(create_md5_file) && select_first([create_md5_file])) then "--CREATE_MD5_FILE" else ""} \
         ~{if (defined(ga4gh_client_secrets) && select_first([ga4gh_client_secrets])) then "--GA4GH_CLIENT_SECRETS" else ""} \
         ~{if (defined(help) && select_first([help])) then "--help" else ""} \
         ~{if defined(max_records_in_ram) then ("--MAX_RECORDS_IN_RAM " + max_records_in_ram) else ''} \
         ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
         ~{if defined(tmp_dir) then ("--TMP_DIR '" + tmp_dir + "'") else ""} \
         ~{if (defined(use_jdk_deflater) && select_first([use_jdk_deflater])) then "--USE_JDK_DEFLATER" else ""} \
         ~{if (defined(use_jdk_inflater) && select_first([use_jdk_inflater])) then "--USE_JDK_INFLATER" else ""} \
         ~{if (defined(validation_stringency) && select_first([validation_stringency])) then "--VALIDATION_STRINGENCY" else ""} \
         ~{if (defined(verbosity) && select_first([verbosity])) then "--VERBOSITY" else ""} \
         ~{if (defined(version) && select_first([version])) then "--version" else ""} \
         ~{if (defined(showhidden) && select_first([showhidden])) then "--showHidden" else ""}
       if [ -f $(echo '~{select_first([outputFilename, "~{basename(inp, ".bam")}.bam"])}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{select_first([outputFilename, "~{basename(inp, ".bam")}.bam"])}' | sed 's/\.[^.]*$//').bai $(echo '~{select_first([outputFilename, "~{basename(inp, ".bam")}.bam"])}' ).bai; fi
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.4.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "~{basename(inp, ".bam")}.bam"])
       File out_bai = select_first([outputFilename, "~{basename(inp, ".bam")}.bam"]) + ".bai"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'Gatk4: ReorderSam'
   doc: |2-

     USAGE: ReorderSam [arguments]

     Not to be confused with SortSam which sorts a SAM or BAM file with a valid sequence dictionary, 
     ReorderSam reorders
     reads in a SAM/BAM file to match the contig ordering in a provided reference file, 
     as determined by exact name matching
     of contigs.  Reads mapped to contigs absent in the new reference 
     are dropped. Runs substantially faster if the input is
     an indexed BAM file.

     Example:

     .. code-tool: none

        java -jar picard.jar ReorderSam \
            INPUT=sample.bam \
            OUTPUT=reordered.bam \
            REFERENCE=reference_with_different_order.fasta
         
     Version:4.1.3.0

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
   - id: inp
     label: inp
     doc: (-I) Input file (SAM or BAM) to extract reads from. Required.
     type: File
     inputBinding:
       prefix: --INPUT
       separate: true
   - id: outputFilename
     label: outputFilename
     doc: (-O) Output file (SAM or BAM) to write extracted reads to. Required.
     type:
     - string
     - 'null'
     default: generated.bam
     inputBinding:
       prefix: --OUTPUT
       valueFrom: $(inputs.inp.basename.replace(/.bam$/, "")).bam
       separate: true
   - id: sequence_dictionary
     label: sequence_dictionary
     doc: |-
       A Sequence Dictionary for the OUTPUT file (can be read from one of the following file types (SAM, BAM, VCF, BCF, Interval List, Fasta, or Dict)
     type: File
     inputBinding:
       prefix: --SEQUENCE_DICTIONARY
       separate: true
   - id: allow_contig_length_discordance
     label: allow_contig_length_discordance
     doc: |-
       (-U)  If true, then permits mapping from a read contig to a new reference contig with the same name but a different length.  Highly dangerous, only use if you know what you are doing.  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ALLOW_CONTIG_LENGTH_DISCORDANCE
       separate: true
   - id: allow_incomplete_dict_concordance
     label: allow_incomplete_dict_concordance
     doc: |-
       (-S)  If true, then allows only a partial overlap of the original contigs with the new reference sequence contigs.  By default, this tool requires a corresponding contig in the new reference for each read contig  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ALLOW_INCOMPLETE_DICT_CONCORDANCE
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
     doc: BAM to write extracted reads to
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
       glob: $(inputs.inp.basename.replace(/.bam$/, "")).bam
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - ReorderSam
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: GatkReorderSam


