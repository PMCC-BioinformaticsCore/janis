:orphan:

Gatk4: AddOrReplaceReadGroups
==========================================================

``GatkAddOrReplaceReadGroups`` · *1 contributor · 4 versions*

USAGE: AddOrReplaceReadGroups [arguments]"
Assigns all the reads in a file to a single new read-group.
This tool accepts INPUT BAM and SAM files or URLs from the <a href='http://ga4gh.org/#/documentation'>Global Alliance
for Genomics and Health (GA4GH)</a>.

Usage example:
++++++++++++++++

.. code-tool: none
   
   java -jar picard.jar AddOrReplaceReadGroups \
      I=input.bam \
      O=output.bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=ILLUMINA \
      RGPU=unit1 \
      RGSM=20
      
Caveats
+++++++++

The value of the tags must adhere (according to the 
<ahref='https://samtools.github.io/hts-specs/SAMv1.pdf'>SAM-spec</a>) with the regex 
<code>'^[ -~]+$'</code> (one or more
characters from the ASCII range 32 through 126). 
In particular &lt;Space&gt; is the only non-printing character allowed.
The program enables only the wholesale assignment of all the reads in the INPUT to a 
single read-group. If your file
already has reads assigned to multiple read-groups, 
the original RG value will be lost. 
For more information about read-groups, see the 
<a href='https://www.broadinstitute.org/gatk/guide/article?id=6472'>GATK Dictionary entry.</a>

Version:4.1.3.0


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.addorreplacereadgroups.versions import Gatk4AddOrReplaceReadGroups_4_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatkaddorreplacereadgroups_step",
           Gatk4AddOrReplaceReadGroups_4_1_2(
               inp=None,
               rglb=None,
               rgpl=None,
               rgpu=None,
               rgsm=None,
           )
       )
       wf.output("out", source=gatkaddorreplacereadgroups_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for GatkAddOrReplaceReadGroups:

.. code-block:: bash

   # user inputs
   janis inputs GatkAddOrReplaceReadGroups > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inp: inp.bam
       rglb: <value>
       rgpl: <value>
       rgpu: <value>
       rgsm: <value>




5. Run GatkAddOrReplaceReadGroups with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       GatkAddOrReplaceReadGroups





Information
------------

:ID: ``GatkAddOrReplaceReadGroups``
:URL: *No URL to the documentation was provided*
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.2.0
:Authors: illusional
:Citations: None
:Created: 2020-05-15
:Updated: 2020-05-15


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
inp                    BAM                      --INPUT                              (-I) Input file (BAM or SAM or a GA4GH url). Required.
rglb                   String                   --RGLB                               (-LB) Read-Group library Required.
rgpl                   String                   --RGPL                               (-PL) Read-Group platform (e.g. ILLUMINA, SOLID) Required.
rgpu                   String                   --RGPU                               (-PU) Read-Group platform unit (eg. run barcode) Required.
rgsm                   String                   --RGSM                               (-SM) Read-Group sample name Required.
javaOptions            Optional<Array<String>>
compression_level      Optional<Integer>                                             Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
outputFilename         Optional<Filename>       --OUTPUT                             (-O) Output file (BAM or SAM). Required.
arguments_file         Optional<File>           --arguments_file                     read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null.
create_index           Optional<Boolean>        --CREATE_INDEX                       Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value: false. Possible values: {true, false}
create_md5_file        Optional<Boolean>        --CREATE_MD5_FILE                    Whether to create an MD5 digest for any BAM or FASTQ files created. Default value: false. Possible values: {true, false}
ga4gh_client_secrets   Optional<Boolean>        --GA4GH_CLIENT_SECRETS               Default value: client_secrets.json.
help                   Optional<Boolean>        --help                               (-h) display the help message Default value: false. Possible values: {true, false}
max_records_in_ram     Optional<Integer>        --MAX_RECORDS_IN_RAM                 When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.  Default value: 500000.
quiet                  Optional<Boolean>        --QUIET                              Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false}
reference_sequence     Optional<File>           --REFERENCE_SEQUENCE                 (-R) Reference sequence file. Default value: null.
rgcn                   Optional<String>         --RGCN                               (-CN) Read-Group sequencing center name Default value: null.
rgds                   Optional<String>         --RGDS                               (-DS) Read-Group description Default value: null.
rgdt                   Optional<Boolean>        --RGDT                               (-DT) Read-Group run date Default value: null.
rgfo                   Optional<String>         --RGFO                               (-FO) Read-Group flow order Default value: null.
rgid                   Optional<String>         --RGID                               (-ID) Read-Group ID Default value: 1.
rgks                   Optional<String>         --RGKS                               (-KS) Read-Group key sequence Default value: null.
rgpg                   Optional<String>         --RGPG                               (-PG) Read-Group program group Default value: null.
rgpi                   Optional<Integer>        --RGPI                               (-PI) Read-Group predicted insert size Default value: null.
rgpm                   Optional<String>         --RGPM                               (-PM) Read-Group platform model Default value: null.
sort_order             Optional<String>         --SORT_ORDER                         (-SO) Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT. Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate, unknown}
tmp_dir                Optional<File>           --TMP_DIR                            One or more directories with space available to be used by this program for temporary storage of working files  This argument may be specified 0 or more times. Default value: null.
use_jdk_deflater       Optional<Boolean>        --USE_JDK_DEFLATER                   (-use_jdk_deflater)  Use the JDK Deflater instead of the Intel Deflater for writing compressed output  Default value: false. Possible values: {true, false}
use_jdk_inflater       Optional<Boolean>        --USE_JDK_INFLATER                   (-use_jdk_inflater)  Use the JDK Inflater instead of the Intel Inflater for reading compressed input  Default value: false. Possible values: {true, false}
validation_stringency  Optional<String>         --VALIDATION_STRINGENCY              Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. Possible values: {STRICT, LENIENT, SILENT}
verbosity              Optional<Boolean>        --VERBOSITY                          Control verbosity of logging. Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG}
version                Optional<Boolean>        --version                            display the version number for this tool Default value: false. Possible values: {true, false}
showhidden             Optional<Boolean>        --showHidden                         (-showHidden)  display hidden arguments  Default value: false. Possible values: {true, false}
=====================  =======================  =======================  ==========  ============================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task GatkAddOrReplaceReadGroups {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       File inp
       String? outputFilename
       String rglb
       String rgpl
       String rgpu
       String rgsm
       File? arguments_file
       Boolean? create_index
       Boolean? create_md5_file
       Boolean? ga4gh_client_secrets
       Boolean? help
       Int? max_records_in_ram
       Boolean? quiet
       File? reference_sequence
       String? rgcn
       String? rgds
       Boolean? rgdt
       String? rgfo
       String? rgid
       String? rgks
       String? rgpg
       Int? rgpi
       String? rgpm
       String? sort_order
       File? tmp_dir
       Boolean? use_jdk_deflater
       Boolean? use_jdk_inflater
       String? validation_stringency
       Boolean? verbosity
       Boolean? version
       Boolean? showhidden
     }
     command <<<
       set -e
       gatk AddOrReplaceReadGroups \
         --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         --INPUT '~{inp}' \
         --OUTPUT '~{select_first([outputFilename, "~{basename(inp, ".bam")}.bam"])}' \
         --RGLB '~{rglb}' \
         --RGPL '~{rgpl}' \
         --RGPU '~{rgpu}' \
         --RGSM '~{rgsm}' \
         ~{if defined(arguments_file) then ("--arguments_file '" + arguments_file + "'") else ""} \
         ~{if (defined(create_index) && select_first([create_index])) then "--CREATE_INDEX" else ""} \
         ~{if (defined(create_md5_file) && select_first([create_md5_file])) then "--CREATE_MD5_FILE" else ""} \
         ~{if (defined(ga4gh_client_secrets) && select_first([ga4gh_client_secrets])) then "--GA4GH_CLIENT_SECRETS" else ""} \
         ~{if (defined(help) && select_first([help])) then "--help" else ""} \
         ~{if defined(max_records_in_ram) then ("--MAX_RECORDS_IN_RAM " + max_records_in_ram) else ''} \
         ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
         ~{if defined(reference_sequence) then ("--REFERENCE_SEQUENCE '" + reference_sequence + "'") else ""} \
         ~{if defined(rgcn) then ("--RGCN '" + rgcn + "'") else ""} \
         ~{if defined(rgds) then ("--RGDS '" + rgds + "'") else ""} \
         ~{if (defined(rgdt) && select_first([rgdt])) then "--RGDT" else ""} \
         ~{if defined(rgfo) then ("--RGFO '" + rgfo + "'") else ""} \
         ~{if defined(rgid) then ("--RGID '" + rgid + "'") else ""} \
         ~{if defined(rgks) then ("--RGKS '" + rgks + "'") else ""} \
         ~{if defined(rgpg) then ("--RGPG '" + rgpg + "'") else ""} \
         ~{if defined(rgpi) then ("--RGPI " + rgpi) else ''} \
         ~{if defined(rgpm) then ("--RGPM '" + rgpm + "'") else ""} \
         ~{if defined(sort_order) then ("--SORT_ORDER '" + sort_order + "'") else ""} \
         ~{if defined(tmp_dir) then ("--TMP_DIR '" + tmp_dir + "'") else ""} \
         ~{if (defined(use_jdk_deflater) && select_first([use_jdk_deflater])) then "--USE_JDK_DEFLATER" else ""} \
         ~{if (defined(use_jdk_inflater) && select_first([use_jdk_inflater])) then "--USE_JDK_INFLATER" else ""} \
         ~{if defined(validation_stringency) then ("--VALIDATION_STRINGENCY '" + validation_stringency + "'") else ""} \
         ~{if (defined(verbosity) && select_first([verbosity])) then "--VERBOSITY" else ""} \
         ~{if (defined(version) && select_first([version])) then "--version" else ""} \
         ~{if (defined(showhidden) && select_first([showhidden])) then "--showHidden" else ""}
       if [ -f $(echo '~{select_first([outputFilename, "~{basename(inp, ".bam")}.bam"])}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{select_first([outputFilename, "~{basename(inp, ".bam")}.bam"])}' | sed 's/\.[^.]*$//').bai $(echo '~{select_first([outputFilename, "~{basename(inp, ".bam")}.bam"])}' ).bai; fi
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
   label: 'Gatk4: AddOrReplaceReadGroups'
   doc: |-
     USAGE: AddOrReplaceReadGroups [arguments]"
     Assigns all the reads in a file to a single new read-group.
     This tool accepts INPUT BAM and SAM files or URLs from the <a href='http://ga4gh.org/#/documentation'>Global Alliance
     for Genomics and Health (GA4GH)</a>.

     Usage example:
     ++++++++++++++++

     .. code-tool: none
     
        java -jar picard.jar AddOrReplaceReadGroups \
           I=input.bam \
           O=output.bam \
           RGID=4 \
           RGLB=lib1 \
           RGPL=ILLUMINA \
           RGPU=unit1 \
           RGSM=20
        
     Caveats
     +++++++++

     The value of the tags must adhere (according to the 
     <ahref='https://samtools.github.io/hts-specs/SAMv1.pdf'>SAM-spec</a>) with the regex 
     <code>'^[ -~]+$'</code> (one or more
     characters from the ASCII range 32 through 126). 
     In particular &lt;Space&gt; is the only non-printing character allowed.
     The program enables only the wholesale assignment of all the reads in the INPUT to a 
     single read-group. If your file
     already has reads assigned to multiple read-groups, 
     the original RG value will be lost. 
     For more information about read-groups, see the 
     <a href='https://www.broadinstitute.org/gatk/guide/article?id=6472'>GATK Dictionary entry.</a>

     Version:4.1.3.0

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
   - id: inp
     label: inp
     doc: (-I) Input file (BAM or SAM or a GA4GH url). Required.
     type: File
     inputBinding:
       prefix: --INPUT
       separate: true
   - id: outputFilename
     label: outputFilename
     doc: (-O) Output file (BAM or SAM). Required.
     type:
     - string
     - 'null'
     default: generated.bam
     inputBinding:
       prefix: --OUTPUT
       valueFrom: $(inputs.inp.basename.replace(/.bam$/, "")).bam
       separate: true
   - id: rglb
     label: rglb
     doc: (-LB) Read-Group library Required.
     type: string
     inputBinding:
       prefix: --RGLB
       separate: true
   - id: rgpl
     label: rgpl
     doc: (-PL) Read-Group platform (e.g. ILLUMINA, SOLID) Required.
     type: string
     inputBinding:
       prefix: --RGPL
       separate: true
   - id: rgpu
     label: rgpu
     doc: (-PU) Read-Group platform unit (eg. run barcode) Required.
     type: string
     inputBinding:
       prefix: --RGPU
       separate: true
   - id: rgsm
     label: rgsm
     doc: (-SM) Read-Group sample name Required.
     type: string
     inputBinding:
       prefix: --RGSM
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
     type:
     - boolean
     - 'null'
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
   - id: rgcn
     label: rgcn
     doc: '(-CN) Read-Group sequencing center name Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --RGCN
       separate: true
   - id: rgds
     label: rgds
     doc: '(-DS) Read-Group description Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --RGDS
       separate: true
   - id: rgdt
     label: rgdt
     doc: '(-DT) Read-Group run date Default value: null.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --RGDT
       separate: true
   - id: rgfo
     label: rgfo
     doc: '(-FO) Read-Group flow order Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --RGFO
       separate: true
   - id: rgid
     label: rgid
     doc: '(-ID) Read-Group ID Default value: 1.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --RGID
       separate: true
   - id: rgks
     label: rgks
     doc: '(-KS) Read-Group key sequence Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --RGKS
       separate: true
   - id: rgpg
     label: rgpg
     doc: '(-PG) Read-Group program group Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --RGPG
       separate: true
   - id: rgpi
     label: rgpi
     doc: '(-PI) Read-Group predicted insert size Default value: null.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --RGPI
       separate: true
   - id: rgpm
     label: rgpm
     doc: '(-PM) Read-Group platform model Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --RGPM
       separate: true
   - id: sort_order
     label: sort_order
     doc: |-
       (-SO) Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT. Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate, unknown} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --SORT_ORDER
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
     - string
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
       glob: $(inputs.inp.basename.replace(/.bam$/, "")).bam
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - AddOrReplaceReadGroups
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: GatkAddOrReplaceReadGroups


