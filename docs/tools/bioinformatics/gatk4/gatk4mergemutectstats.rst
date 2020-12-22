:orphan:

GATK4: MergeMutectStats
===============================================

``Gatk4MergeMutectStats`` · *1 contributor · 6 versions*

TBD


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.mergemutectstats.versions import Gatk4MergeMutectStats_4_1_8

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4mergemutectstats_step",
           Gatk4MergeMutectStats_4_1_8(
               statsFiles=None,
           )
       )
       wf.output("out", source=gatk4mergemutectstats_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4MergeMutectStats:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4MergeMutectStats > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       statsFiles:
       - statsFiles_0.txt
       - statsFiles_1.txt




5. Run Gatk4MergeMutectStats with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4MergeMutectStats





Information
------------

:ID: ``Gatk4MergeMutectStats``
:URL: `TBD <TBD>`_
:Versions: 4.1.8.1, 4.1.7.0, 4.1.6.0, 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.8.1
:Authors: Hollizeck Sebastian
:Citations: TBD
:Created: 2019-09-09
:Updated: 2019-09-09


Outputs
-----------

======  ========  ========================
name    type      documentation
======  ========  ========================
out     TextFile  Merged callability stats
======  ========  ========================


Additional configuration (inputs)
---------------------------------

=================  =======================  ========  ==========  ========================================================================================
name               type                     prefix      position  documentation
=================  =======================  ========  ==========  ========================================================================================
statsFiles         Array<TextFile>          --stats            0  Callability stats
javaOptions        Optional<Array<String>>
compression_level  Optional<Integer>                              Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
mergedStatsOut     Optional<Filename>       -O                 1
=================  =======================  ========  ==========  ========================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4MergeMutectStats {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       Array[File] statsFiles
       String? mergedStatsOut
     }
     command <<<
       set -e
       gatk MergeMutectStats \
         --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         ~{if length(statsFiles) > 0 then "--stats '" + sep("' --stats '", statsFiles) + "'" else ""} \
         -O '~{select_first([mergedStatsOut, "generated.txt"])}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.8.1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([mergedStatsOut, "generated.txt"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: MergeMutectStats'
   doc: TBD

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.8.1

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
   - id: statsFiles
     label: statsFiles
     doc: Callability stats
     type:
       type: array
       inputBinding:
         prefix: --stats
       items: File
     inputBinding:
       position: 0
   - id: mergedStatsOut
     label: mergedStatsOut
     type:
     - string
     - 'null'
     default: generated.txt
     inputBinding:
       prefix: -O
       position: 1

   outputs:
   - id: out
     label: out
     doc: Merged callability stats
     type: File
     outputBinding:
       glob: generated.txt
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - MergeMutectStats
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4MergeMutectStats


