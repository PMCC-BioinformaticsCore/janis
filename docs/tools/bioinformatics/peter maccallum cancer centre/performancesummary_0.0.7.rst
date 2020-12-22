:orphan:

Performance Summary
========================================

``performanceSummary`` · *1 contributor · 1 version*

usage: performance_summary.py [-h] --flagstat FLAGSTAT
                              --collect_insert_metrics COLLECT_INSERT_METRICS
                              --coverage COVERAGE -o O
                              [--target_flagstat TARGET_FLAGSTAT]
                              [--rmdup_flagstat RMDUP_FLAGSTAT] [--genome]

Performance summary of bam

required arguments:
  --flagstat FLAGSTAT   output of samtools flagstat on bam
  --collect_insert_metrics COLLECT_INSERT_METRICS
                        output of CollectInsertMetrics (GATK or Picard) on bam
  --coverage COVERAGE   output of bedtools coverageBed for targeted bam;
                        bedtools genomeCoverageBed for whole genome bam
  -o O                  output summary csv name

optional arguments:
  -h, --help            show this help message and exit
  --target_flagstat TARGET_FLAGSTAT
                        output of samtools flagstat of bam target on target
                        bed. Only specified for targeted bam
  --rmdup_flagstat RMDUP_FLAGSTAT
                        output of samtools flagstat of removed duplicates bam.
                        File to be used to extract mapping infomation if
                        specified, instead of the --flagstat file.
  --genome              calculate statistics for whole genome data.
                        --target_flagstat must not be speicified
        


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.performancesummary.versions import PerformanceSummary_0_0_7

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "performancesummary_step",
           PerformanceSummary_0_0_7(
               flagstat=None,
               collectInsertSizeMetrics=None,
               coverage=None,
           )
       )
       wf.output("out", source=performancesummary_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for performanceSummary:

.. code-block:: bash

   # user inputs
   janis inputs performanceSummary > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       collectInsertSizeMetrics: collectInsertSizeMetrics
       coverage: coverage
       flagstat: flagstat




5. Run performanceSummary with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       performanceSummary





Information
------------

:ID: ``performanceSummary``
:URL: `https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/performance <https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/performance>`_
:Versions: 0.0.7
:Container: michaelfranklin/pmacutil:0.0.7
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-04-03 00:00:00
:Updated: 2020-04-03 00:00:00


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     csv
======  ======  ===============


Additional configuration (inputs)
---------------------------------

========================  ==================  ========================  ==========  ==================================================================================================================================================
name                      type                prefix                    position    documentation
========================  ==================  ========================  ==========  ==================================================================================================================================================
flagstat                  File                --flagstat                            output of samtools flagstat on bam
collectInsertSizeMetrics  File                --collect_insert_metrics              output of CollectInsertMetrics (GATK or Picard) on bam
coverage                  File                --coverage                            output of bedtools coverageBed for targeted bam; bedtools genomeCoverageBed for whole genome bam
outputPrefix              Optional<Filename>  -o                                    prefix of output summary csv
targetFlagstat            Optional<File>      --target_flagstat                     output of samtools flagstat of bam target on target bed. Only specified for targeted bam
rmdupFlagstat             Optional<File>      --rmdup_flagstat                      output of samtools flagstat of removed duplicates bam. File to be used to extract mapping infomation if specified, instead of the --flagstat file.
genome                    Optional<Boolean>   --genome                              calculate statistics for whole genome data.--target_flagstat must not be speicified
========================  ==================  ========================  ==========  ==================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task performanceSummary {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File flagstat
       File collectInsertSizeMetrics
       File coverage
       String? outputPrefix
       File? targetFlagstat
       File? rmdupFlagstat
       Boolean? genome
     }
     command <<<
       set -e
       performance_summary.py \
         --flagstat '~{flagstat}' \
         --collect_insert_metrics '~{collectInsertSizeMetrics}' \
         --coverage '~{coverage}' \
         -o '~{select_first([outputPrefix, "generated.csv"])}' \
         ~{if defined(targetFlagstat) then ("--target_flagstat '" + targetFlagstat + "'") else ""} \
         ~{if defined(rmdupFlagstat) then ("--rmdup_flagstat '" + rmdupFlagstat + "'") else ""} \
         ~{if (defined(genome) && select_first([genome])) then "--genome" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "michaelfranklin/pmacutil:0.0.7"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = (select_first([outputPrefix, "generated.csv"]) + ".csv")
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Performance Summary
   doc: |-
     usage: performance_summary.py [-h] --flagstat FLAGSTAT
                                   --collect_insert_metrics COLLECT_INSERT_METRICS
                                   --coverage COVERAGE -o O
                                   [--target_flagstat TARGET_FLAGSTAT]
                                   [--rmdup_flagstat RMDUP_FLAGSTAT] [--genome]

     Performance summary of bam

     required arguments:
       --flagstat FLAGSTAT   output of samtools flagstat on bam
       --collect_insert_metrics COLLECT_INSERT_METRICS
                             output of CollectInsertMetrics (GATK or Picard) on bam
       --coverage COVERAGE   output of bedtools coverageBed for targeted bam;
                             bedtools genomeCoverageBed for whole genome bam
       -o O                  output summary csv name

     optional arguments:
       -h, --help            show this help message and exit
       --target_flagstat TARGET_FLAGSTAT
                             output of samtools flagstat of bam target on target
                             bed. Only specified for targeted bam
       --rmdup_flagstat RMDUP_FLAGSTAT
                             output of samtools flagstat of removed duplicates bam.
                             File to be used to extract mapping infomation if
                             specified, instead of the --flagstat file.
       --genome              calculate statistics for whole genome data.
                             --target_flagstat must not be speicified
          

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: michaelfranklin/pmacutil:0.0.7

   inputs:
   - id: flagstat
     label: flagstat
     doc: output of samtools flagstat on bam
     type: File
     inputBinding:
       prefix: --flagstat
   - id: collectInsertSizeMetrics
     label: collectInsertSizeMetrics
     doc: output of CollectInsertMetrics (GATK or Picard) on bam
     type: File
     inputBinding:
       prefix: --collect_insert_metrics
   - id: coverage
     label: coverage
     doc: |-
       output of bedtools coverageBed for targeted bam; bedtools genomeCoverageBed for whole genome bam
     type: File
     inputBinding:
       prefix: --coverage
   - id: outputPrefix
     label: outputPrefix
     doc: prefix of output summary csv
     type:
     - string
     - 'null'
     default: generated.csv
     inputBinding:
       prefix: -o
   - id: targetFlagstat
     label: targetFlagstat
     doc: |-
       output of samtools flagstat of bam target on target bed. Only specified for targeted bam
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --target_flagstat
   - id: rmdupFlagstat
     label: rmdupFlagstat
     doc: |-
       output of samtools flagstat of removed duplicates bam. File to be used to extract mapping infomation if specified, instead of the --flagstat file.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --rmdup_flagstat
   - id: genome
     label: genome
     doc: |-
       calculate statistics for whole genome data.--target_flagstat must not be speicified
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --genome

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".csv"))
       outputEval: $((inputs.outputPrefix.basename + ".csv"))
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: performance_summary.py
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: performanceSummary


