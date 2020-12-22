:orphan:

GATK4: GatherBQSRReports
=================================================

``Gatk4GatherBQSRReports`` · *1 contributor · 4 versions*

USAGE: GatherBQSRReports [arguments]
Gathers scattered BQSR recalibration reports into a single file
Version:4.1.3.0



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.gatherbqsrreports.versions import Gatk4GatherBQSRReports_4_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4gatherbqsrreports_step",
           Gatk4GatherBQSRReports_4_1_2(

           )
       )
       wf.output("out", source=gatk4gatherbqsrreports_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4GatherBQSRReports:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4GatherBQSRReports > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run Gatk4GatherBQSRReports with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4GatherBQSRReports





Information
------------

:ID: ``Gatk4GatherBQSRReports``
:URL: *No URL to the documentation was provided*
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.2.0
:Authors: Michael Franklin
:Citations: None
:Created: 2020-05-18
:Updated: 2020-05-18


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     tsv
======  ======  ===============


Additional configuration (inputs)
---------------------------------

==========================  =======================  ================================  ==========  ======================================================================================================================================
name                        type                     prefix                            position    documentation
==========================  =======================  ================================  ==========  ======================================================================================================================================
javaOptions                 Optional<Array<String>>
compression_level           Optional<Integer>                                                      Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
reports                     Optional<Array<tsv>>     --input                                       (-I) List of scattered BQSR report files This argument must be specified at least once. Required.
outputFilename              Optional<Filename>       --output                                      (-O) File to output the gathered file to Required.
arguments_file              Optional<File>           --arguments_file                              read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null.
gatkConfigFile              Optional<String>         --gatk-config-file                            A configuration file to use with the GATK. Default value: null.
gcsMaxRetries               Optional<Integer>        --gcs-max-retries                             (-gcs-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20.
gcsProjectForRequesterPays  Optional<String>         --gcs-project-for-requester-pays              Project to bill when accessing 'requester pays' buckets. If unset, these buckets cannot be accessed.  Default value: .
help                        Optional<Boolean>        --help                                        (-h) display the help message Default value: false. Possible values: {true, false}
quiet                       Optional<Boolean>        --QUIET                                       Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false}
tmpDir                      Optional<Boolean>        --tmp-dir                                     Temp directory to use. Default value: null.
useJdkDeflater              Optional<Boolean>        --use-jdk-deflater                            (-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false}
useJdkInflater              Optional<Boolean>        --use-jdk-inflater                            (-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false}
verbosity                   Optional<Boolean>        --verbosity                                   (-verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG}
version                     Optional<Boolean>        --version                                     display the version number for this tool Default value: false. Possible values: {true, false}
showhidden                  Optional<Boolean>        --showHidden                                  (-showHidden)  display hidden arguments  Default value: false. Possible values: {true, false}
==========================  =======================  ================================  ==========  ======================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4GatherBQSRReports {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       Array[File]? reports
       String? outputFilename
       File? arguments_file
       String? gatkConfigFile
       Int? gcsMaxRetries
       String? gcsProjectForRequesterPays
       Boolean? help
       Boolean? quiet
       Boolean? tmpDir
       Boolean? useJdkDeflater
       Boolean? useJdkInflater
       Boolean? verbosity
       Boolean? version
       Boolean? showhidden
     }
     command <<<
       set -e
       gatk GatherBQSRReports \
         --java-options '-Xmx~{((select_first([runtime_memory, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         ~{if (defined(reports) && length(select_first([reports])) > 0) then "--input '" + sep("' --input '", select_first([reports])) + "'" else ""} \
         --output '~{select_first([outputFilename, "generated.recal_data.tsv"])}' \
         ~{if defined(arguments_file) then ("--arguments_file '" + arguments_file + "'") else ""} \
         ~{if defined(gatkConfigFile) then ("--gatk-config-file '" + gatkConfigFile + "'") else ""} \
         ~{if defined(gcsMaxRetries) then ("--gcs-max-retries " + gcsMaxRetries) else ''} \
         ~{if defined(gcsProjectForRequesterPays) then ("--gcs-project-for-requester-pays '" + gcsProjectForRequesterPays + "'") else ""} \
         ~{if (defined(help) && select_first([help])) then "--help" else ""} \
         ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
         ~{if (defined(tmpDir) && select_first([tmpDir])) then "--tmp-dir" else ""} \
         ~{if (defined(useJdkDeflater) && select_first([useJdkDeflater])) then "--use-jdk-deflater" else ""} \
         ~{if (defined(useJdkInflater) && select_first([useJdkInflater])) then "--use-jdk-inflater" else ""} \
         ~{if (defined(verbosity) && select_first([verbosity])) then "--verbosity" else ""} \
         ~{if (defined(version) && select_first([version])) then "--version" else ""} \
         ~{if (defined(showhidden) && select_first([showhidden])) then "--showHidden" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.2.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.recal_data.tsv"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: GatherBQSRReports'
   doc: |
     USAGE: GatherBQSRReports [arguments]
     Gathers scattered BQSR recalibration reports into a single file
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
   - id: reports
     label: reports
     doc: |-
       (-I) List of scattered BQSR report files This argument must be specified at least once. Required. 
     type:
     - type: array
       inputBinding:
         prefix: --input
         separate: true
       items: File
     - 'null'
     inputBinding: {}
   - id: outputFilename
     label: outputFilename
     doc: (-O) File to output the gathered file to Required.
     type:
     - string
     - 'null'
     default: generated.recal_data.tsv
     inputBinding:
       prefix: --output
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
   - id: gatkConfigFile
     label: gatkConfigFile
     doc: 'A configuration file to use with the GATK. Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --gatk-config-file
       separate: true
   - id: gcsMaxRetries
     label: gcsMaxRetries
     doc: |-
       (-gcs-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --gcs-max-retries
       separate: true
   - id: gcsProjectForRequesterPays
     label: gcsProjectForRequesterPays
     doc: |2-
        Project to bill when accessing 'requester pays' buckets. If unset, these buckets cannot be accessed.  Default value: . 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --gcs-project-for-requester-pays
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
   - id: tmpDir
     label: tmpDir
     doc: 'Temp directory to use. Default value: null.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --tmp-dir
       separate: true
   - id: useJdkDeflater
     label: useJdkDeflater
     doc: |-
       (-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --use-jdk-deflater
       separate: true
   - id: useJdkInflater
     label: useJdkInflater
     doc: |-
       (-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --use-jdk-inflater
       separate: true
   - id: verbosity
     label: verbosity
     doc: |-
       (-verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --verbosity
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
     outputBinding:
       glob: generated.recal_data.tsv
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - GatherBQSRReports
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4GatherBQSRReports


