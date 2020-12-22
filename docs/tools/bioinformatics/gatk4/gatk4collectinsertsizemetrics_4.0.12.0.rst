:orphan:

GATK4: CollectInsertSizeMetrics
===============================================================

``Gatk4CollectInsertSizeMetrics`` · *1 contributor · 4 versions*

Provides useful metrics for validating library construction including the insert size distribution and read orientation of paired-end libraries


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.collectinsertsizemetrics.versions import Gatk4CollectInsertSizeMetrics_4_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4collectinsertsizemetrics_step",
           Gatk4CollectInsertSizeMetrics_4_0(
               bam=None,
           )
       )
       wf.output("out", source=gatk4collectinsertsizemetrics_step.out)
       wf.output("outHistogram", source=gatk4collectinsertsizemetrics_step.outHistogram)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4CollectInsertSizeMetrics:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4CollectInsertSizeMetrics > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam




5. Run Gatk4CollectInsertSizeMetrics with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4CollectInsertSizeMetrics





Information
------------

:ID: ``Gatk4CollectInsertSizeMetrics``
:URL: `https://gatk.broadinstitute.org/hc/en-us/articles/360036715591-CollectInsertSizeMetrics-Picard- <https://gatk.broadinstitute.org/hc/en-us/articles/360036715591-CollectInsertSizeMetrics-Picard->`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.0.12.0
:Authors: Jiaan Yu
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2020-02-17
:Updated: 2020-02-17


Outputs
-----------

============  ========  ===============
name          type      documentation
============  ========  ===============
out           TextFile
outHistogram  File
============  ========  ===============


Additional configuration (inputs)
---------------------------------

=======================  =======================  ===========================  ==========  ============================================================================================================================================================================================================================================================================================================================
name                     type                     prefix                         position  documentation
=======================  =======================  ===========================  ==========  ============================================================================================================================================================================================================================================================================================================================
bam                      IndexedBam               -I                                   10  Input SAM or BAM file.  Required.
javaOptions              Optional<Array<String>>
compression_level        Optional<Integer>                                                 Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
outputFilename           Optional<Filename>       -O                                       File to write the output to.  Required.
outputHistogram          Optional<Filename>       -H                                       File to write insert size Histogram chart to.  Required.
argumentsFile            Optional<Array<File>>    --arguments_file                     10  read one or more arguments files and add them to the command line
assumeSorted             Optional<Boolean>        --ASSUME_SORTED                      11  If true (default), then the sort order in the header file will be ignored.  Default value: true. Possible values: {true, false}
deviations               Optional<Double>         --DEVIATIONS                         11  Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. This is done because insert size data typically includes enough anomalous values from chimeras and other artifacts to make the mean and sd grossly misleading regarding the real distribution.  Default value: 10.0.
histogramWidth           Optional<Integer>        --HISTOGRAM_WIDTH                    11  Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. Also, when calculating mean and standard deviation, only bins <= Histogram_WIDTH will be included.  Default value: null.
includeDuplicates        Optional<Boolean>        --INCLUDE_DUPLICATES                 11  If true, also include reads marked as duplicates in the insert size histogram.  Default value: false. Possible values: {true, false}
metricAccumulationLevel  Optional<String>         --METRIC_ACCUMULATION_LEVEL          11  The level(s) at  which to accumulate metrics.    This argument may be specified 0 or more times. Default value: [ALL_READS]. Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} .
minimumPCT               Optional<Float>          --MINIMUM_PCT                        11  When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this percentage of overall reads. (Range: 0 to 1).  Default value: 0.05.
stopAfter                Optional<Integer>        --STOP_AFTER                         11  Stop after  processing N reads, mainly for debugging.  Default value: 0.
version                  Optional<Boolean>        --version                            11  display the version number for this tool Default value: false. Possible values: {true, false}
showHidden               Optional<Boolean>        --showHidden                         11  display hidden  arguments  Default  value: false.  Possible values: {true, false}
=======================  =======================  ===========================  ==========  ============================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4CollectInsertSizeMetrics {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       File bam
       File bam_bai
       String? outputFilename
       String? outputHistogram
       Array[File]? argumentsFile
       Boolean? assumeSorted
       Float? deviations
       Int? histogramWidth
       Boolean? includeDuplicates
       String? metricAccumulationLevel
       Float? minimumPCT
       Int? stopAfter
       Boolean? version
       Boolean? showHidden
     }
     command <<<
       set -e
       gatk CollectInsertSizeMetrics \
         --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         -O '~{select_first([outputFilename, "~{basename(bam, ".bam")}.metrics.txt"])}' \
         -H '~{select_first([outputHistogram, "~{basename(bam, ".bam")}.histogram.pdf"])}' \
         -I '~{bam}' \
         ~{if (defined(argumentsFile) && length(select_first([argumentsFile])) > 0) then "--arguments_file '" + sep("' --arguments_file '", select_first([argumentsFile])) + "'" else ""} \
         ~{if (defined(assumeSorted) && select_first([assumeSorted])) then "--ASSUME_SORTED" else ""} \
         ~{if defined(deviations) then ("--DEVIATIONS " + deviations) else ''} \
         ~{if defined(histogramWidth) then ("--HISTOGRAM_WIDTH " + histogramWidth) else ''} \
         ~{if (defined(includeDuplicates) && select_first([includeDuplicates])) then "--INCLUDE_DUPLICATES" else ""} \
         ~{if defined(metricAccumulationLevel) then ("--METRIC_ACCUMULATION_LEVEL '" + metricAccumulationLevel + "'") else ""} \
         ~{if defined(minimumPCT) then ("--MINIMUM_PCT " + minimumPCT) else ''} \
         ~{if defined(stopAfter) then ("--STOP_AFTER " + stopAfter) else ''} \
         ~{if (defined(version) && select_first([version])) then "--version" else ""} \
         ~{if (defined(showHidden) && select_first([showHidden])) then "--showHidden" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.0.12.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "~{basename(bam, ".bam")}.metrics.txt"])
       File outHistogram = select_first([outputHistogram, "~{basename(bam, ".bam")}.histogram.pdf"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: CollectInsertSizeMetrics'
   doc: |-
     Provides useful metrics for validating library construction including the insert size distribution and read orientation of paired-end libraries

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.0.12.0

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
     doc: Input SAM or BAM file.  Required.
     type: File
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       prefix: -I
       position: 10
   - id: outputFilename
     label: outputFilename
     doc: File to write the output to.  Required.
     type:
     - string
     - 'null'
     default: generated.metrics.txt
     inputBinding:
       prefix: -O
       valueFrom: $(inputs.bam.basename.replace(/.bam$/, "")).metrics.txt
   - id: outputHistogram
     label: outputHistogram
     doc: 'File to write insert size Histogram chart to.  Required. '
     type:
     - string
     - 'null'
     default: generated.histogram.pdf
     inputBinding:
       prefix: -H
       valueFrom: $(inputs.bam.basename.replace(/.bam$/, "")).histogram.pdf
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
   - id: assumeSorted
     label: assumeSorted
     doc: |-
       If true (default), then the sort order in the header file will be ignored.  Default value: true. Possible values: {true, false}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ASSUME_SORTED
       position: 11
   - id: deviations
     label: deviations
     doc: |-
       Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. This is done because insert size data typically includes enough anomalous values from chimeras and other artifacts to make the mean and sd grossly misleading regarding the real distribution.  Default value: 10.0. 
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --DEVIATIONS
       position: 11
   - id: histogramWidth
     label: histogramWidth
     doc: |-
       Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. Also, when calculating mean and standard deviation, only bins <= Histogram_WIDTH will be included.  Default value: null. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --HISTOGRAM_WIDTH
       position: 11
   - id: includeDuplicates
     label: includeDuplicates
     doc: |-
       If true, also include reads marked as duplicates in the insert size histogram.  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --INCLUDE_DUPLICATES
       position: 11
   - id: metricAccumulationLevel
     label: metricAccumulationLevel
     doc: |-
       The level(s) at  which to accumulate metrics.    This argument may be specified 0 or more times. Default value: [ALL_READS]. Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} .
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --METRIC_ACCUMULATION_LEVEL
       position: 11
   - id: minimumPCT
     label: minimumPCT
     doc: |-
       When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this percentage of overall reads. (Range: 0 to 1).  Default value: 0.05.
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --MINIMUM_PCT
       position: 11
   - id: stopAfter
     label: stopAfter
     doc: 'Stop after  processing N reads, mainly for debugging.  Default value: 0. '
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --STOP_AFTER
       position: 11
   - id: version
     label: version
     doc: |-
       display the version number for this tool Default value: false. Possible values: {true, false}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --version
       position: 11
   - id: showHidden
     label: showHidden
     doc: |-
       display hidden  arguments  Default  value: false.  Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --showHidden
       position: 11

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: $(inputs.bam.basename.replace(/.bam$/, "")).metrics.txt
       loadContents: false
   - id: outHistogram
     label: outHistogram
     type: File
     outputBinding:
       glob: $(inputs.bam.basename.replace(/.bam$/, "")).histogram.pdf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - CollectInsertSizeMetrics
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4CollectInsertSizeMetrics


