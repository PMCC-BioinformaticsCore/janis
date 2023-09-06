:orphan:

Refilter Strelka2 Variant Calls
=======================================================

``refilterStrelka2Calls`` · *1 contributor · 1 version*

Usage: filterStrelkaCalls.R [options]



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.dawson.refilterstrelka2calls.versions import RefilterStrelka2Calls_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "refilterstrelka2calls_step",
           RefilterStrelka2Calls_0_1(
               inputFiles=None,
               MQ=None,
               DP=None,
               EVS=None,
               RPRS=None,
               minAD=None,
               threads=None,
               outputFolder=None,
           )
       )
       wf.output("out", source=refilterstrelka2calls_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for refilterStrelka2Calls:

.. code-block:: bash

   # user inputs
   janis inputs refilterStrelka2Calls > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inputFiles:
       - inputFiles_0.vcf.gz
       - inputFiles_1.vcf.gz




5. Run refilterStrelka2Calls with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       refilterStrelka2Calls

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          refilterStrelka2Calls

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``refilterStrelka2Calls``
:URL: *No URL to the documentation was provided*
:Versions: 0.2
:Container: shollizeck/dawsontoolkit:0.2
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2019-10-19
:Updated: 2019-10-25


Outputs
-----------

======  ==========  =================
name    type        documentation
======  ==========  =================
out     Array<VCF>  To determine type
======  ==========  =================


Additional configuration (inputs)
---------------------------------

============  =======================  =============  ==========  ========================================================================
name          type                     prefix         position    documentation
============  =======================  =============  ==========  ========================================================================
inputFiles    Array<Gzipped<VCF>>      -i                         comma seperated list of vcfs
MQ            Integer                  --mq                       minimum mapping quality for a variant to be accepted (default: 15)
DP            Integer                  --dp                       minimum depth of coverage for a variant to be accepted (default: 10)
EVS           Integer                  --evs                      minimum phred scaled evidence for a variant to be accepted (default: 20)
RPRS          Integer                  --rprs                     minimum phred scaled evidence for a variant to be accepted (default: 20)
minAD         Integer                  --minAD                    minimum allelic depth for a variant to be accepted (default: 2)
threads       Integer                  -t                         amount of threads to use for parallelization (default: 5)
outputFolder  String                   -o                         Name of the normal sample (default: infered from all sample names)
interval      Optional<String>         -L                         interval to call on (default: everything)
normalName    Optional<String>         -n                         Name of the normal sample (default: infered from all sample names)
sampleNames   Optional<Array<String>>  --sampleNames              Name of the normal sample (default: infered from all sample names)
============  =======================  =============  ==========  ========================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task refilterStrelka2Calls {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       Array[File] inputFiles
       Array[File] inputFiles_tbi
       Int? MQ
       Int? DP
       Int? EVS
       Int? RPRS
       Int? minAD
       Int? threads
       String? interval
       String? normalName
       Array[String]? sampleNames
       String? outputFolder
     }

     command <<<
       set -e
       filterStrelkaCalls.R \
         ~{if length(inputFiles) > 0 then "-i '" + sep("','", inputFiles) + "'" else ""} \
         --mq ~{select_first([MQ, 15])} \
         --dp ~{select_first([DP, 10])} \
         --evs ~{select_first([EVS, 20])} \
         --rprs ~{select_first([RPRS, -10])} \
         --minAD ~{select_first([minAD, 2])} \
         -t ~{select_first([threads, select_first([runtime_cpu, 1])])} \
         ~{if defined(interval) then ("-L '" + interval + "'") else ""} \
         ~{if defined(normalName) then ("-n '" + normalName + "'") else ""} \
         ~{if (defined(sampleNames) && length(select_first([sampleNames])) > 0) then "--sampleNames '" + sep("','", select_first([sampleNames])) + "'" else ""} \
         -o '~{select_first([outputFolder, "./"])}'
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 20, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "shollizeck/dawsontoolkit:0.2"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 48, 4])}G"
       preemptible: 2
     }

     output {
       Array[File] out = glob("*.refiltered.vcf")
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Refilter Strelka2 Variant Calls

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: shollizeck/dawsontoolkit:0.2

   inputs:
   - id: inputFiles
     label: inputFiles
     doc: comma seperated list of vcfs
     type:
       type: array
       items: File
     inputBinding:
       prefix: -i
       itemSeparator: ','
   - id: MQ
     label: MQ
     doc: 'minimum mapping quality for a variant to be accepted (default: 15)'
     type: int
     default: 15
     inputBinding:
       prefix: --mq
   - id: DP
     label: DP
     doc: 'minimum depth of coverage for a variant to be accepted (default: 10)'
     type: int
     default: 10
     inputBinding:
       prefix: --dp
   - id: EVS
     label: EVS
     doc: 'minimum phred scaled evidence for a variant to be accepted (default: 20)'
     type: int
     default: 20
     inputBinding:
       prefix: --evs
   - id: RPRS
     label: RPRS
     doc: 'minimum phred scaled evidence for a variant to be accepted (default: 20)'
     type: int
     default: -10
     inputBinding:
       prefix: --rprs
   - id: minAD
     label: minAD
     doc: 'minimum allelic depth for a variant to be accepted (default: 2)'
     type: int
     default: 2
     inputBinding:
       prefix: --minAD
   - id: threads
     label: threads
     doc: 'amount of threads to use for parallelization (default: 5)'
     type: int
     inputBinding:
       prefix: -t
       valueFrom: |-
         $([inputs.runtime_cpu, 20, 1].filter(function (inner) { return inner != null })[0])
   - id: interval
     label: interval
     doc: 'interval to call on (default: everything)'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -L
   - id: normalName
     label: normalName
     doc: 'Name of the normal sample (default: infered from all sample names)'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -n
   - id: sampleNames
     label: sampleNames
     doc: 'Name of the normal sample (default: infered from all sample names)'
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --sampleNames
       itemSeparator: ','
   - id: outputFolder
     label: outputFolder
     doc: 'Name of the normal sample (default: infered from all sample names)'
     type: string
     default: ./
     inputBinding:
       prefix: -o

   outputs:
   - id: out
     label: out
     doc: To determine type
     type:
       type: array
       items: File
     outputBinding:
       glob: '*.refiltered.vcf'
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: filterStrelkaCalls.R
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: refilterStrelka2Calls


