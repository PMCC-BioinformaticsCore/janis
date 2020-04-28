:orphan:

Performance Summary
========================================

*1 contributor · 2 versions*

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
:Versions: dev, 0.0.7
:Container: michaelfranklin/pmacutil:0.0.7
:Authors: Jiaan Yu
:Citations: None
:Created: None
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
