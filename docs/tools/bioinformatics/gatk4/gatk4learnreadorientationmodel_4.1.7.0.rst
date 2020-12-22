:orphan:

GATK4: LearnReadOrientationModel
=================================================================

``Gatk4LearnReadOrientationModel`` · *1 contributor · 6 versions*

TBD


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.learnreadorientationmodel.versions import Gatk4LearnReadOrientationModel_4_1_7

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4learnreadorientationmodel_step",
           Gatk4LearnReadOrientationModel_4_1_7(
               f1r2CountsFiles=None,
           )
       )
       wf.output("out", source=gatk4learnreadorientationmodel_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4LearnReadOrientationModel:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4LearnReadOrientationModel > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       f1r2CountsFiles:
       - f1r2CountsFiles_0.tar.gz
       - f1r2CountsFiles_1.tar.gz




5. Run Gatk4LearnReadOrientationModel with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4LearnReadOrientationModel





Information
------------

:ID: ``Gatk4LearnReadOrientationModel``
:URL: `TBD <TBD>`_
:Versions: 4.1.8.1, 4.1.7.0, 4.1.6.0, 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.7.0
:Authors: Hollizeck Sebastian
:Citations: TBD
:Created: 2019-09-09
:Updated: 2019-09-09


Outputs
-----------

======  =================  =============================================================
name    type               documentation
======  =================  =============================================================
out     CompressedTarFile  Model file containing information about fragment orientations
======  =================  =============================================================


Additional configuration (inputs)
---------------------------------

=================  ========================  ===================  ==========  ========================================================================================
name               type                      prefix                 position  documentation
=================  ========================  ===================  ==========  ========================================================================================
f1r2CountsFiles    Array<CompressedTarFile>  -I                            0  Counts for the read orientation of fragments
javaOptions        Optional<Array<String>>
compression_level  Optional<Integer>                                          Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
numEmIterations    Optional<Integer>         --num-em-iterations           1  Amount of iterations for the em process before it bails
modelFileOut       Optional<Filename>        -O                            3
=================  ========================  ===================  ==========  ========================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4LearnReadOrientationModel {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       Array[File] f1r2CountsFiles
       Int? numEmIterations
       String? modelFileOut
     }
     command <<<
       set -e
       gatk LearnReadOrientationModel \
         --java-options '-Xmx~{((select_first([runtime_memory, 32, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         ~{if length(f1r2CountsFiles) > 0 then "-I '" + sep("' -I '", f1r2CountsFiles) + "'" else ""} \
         ~{if defined(select_first([numEmIterations, 30])) then ("--num-em-iterations " + select_first([numEmIterations, 30])) else ''} \
         -O '~{select_first([modelFileOut, "generated.tar.gz"])}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.7.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 32, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([modelFileOut, "generated.tar.gz"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: LearnReadOrientationModel'
   doc: TBD

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.7.0

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
   - id: f1r2CountsFiles
     label: f1r2CountsFiles
     doc: Counts for the read orientation of fragments
     type:
       type: array
       inputBinding:
         prefix: -I
       items: File
     inputBinding:
       position: 0
   - id: numEmIterations
     label: numEmIterations
     doc: Amount of iterations for the em process before it bails
     type: int
     default: 30
     inputBinding:
       prefix: --num-em-iterations
       position: 1
   - id: modelFileOut
     label: modelFileOut
     type:
     - string
     - 'null'
     default: generated.tar.gz
     inputBinding:
       prefix: -O
       position: 3

   outputs:
   - id: out
     label: out
     doc: Model file containing information about fragment orientations
     type: File
     outputBinding:
       glob: generated.tar.gz
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - LearnReadOrientationModel
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 32, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4LearnReadOrientationModel


