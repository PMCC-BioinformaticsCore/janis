:orphan:

GATK4: CalculateContamination
===========================================================

``Gatk4CalculateContamination`` · *1 contributor · 6 versions*

Calculates the fraction of reads coming from cross-sample contamination, given results from GetPileupSummaries. The resulting contamination table is used with FilterMutectCalls.

This tool is featured in the Somatic Short Mutation calling Best Practice Workflow. See Tutorial#11136 for a step-by-step description of the workflow and Article#11127 for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the Mutect2 WDL scripts directory.

This tool borrows from ContEst by Cibulskis et al the idea of estimating contamination from ref reads at hom alt sites. However, ContEst uses a probabilistic model that assumes a diploid genotype with no copy number variation and independent contaminating reads. That is, ContEst assumes that each contaminating read is drawn randomly and independently from a different human. This tool uses a simpler estimate of contamination that relaxes these assumptions. In particular, it works in the presence of copy number variations and with an arbitrary number of contaminating samples. In addition, this tool is designed to work well with no matched normal data. However, one can run GetPileupSummaries on a matched normal bam file and input the result to this tool.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.calculatecontaminations.versions import Gatk4CalculateContamination_4_1_7

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4calculatecontamination_step",
           Gatk4CalculateContamination_4_1_7(
               pileupTable=None,
           )
       )
       wf.output("contOut", source=gatk4calculatecontamination_step.contOut)
       wf.output("segOut", source=gatk4calculatecontamination_step.segOut)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4CalculateContamination:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4CalculateContamination > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       pileupTable: pileupTable




5. Run Gatk4CalculateContamination with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4CalculateContamination





Information
------------

:ID: ``Gatk4CalculateContamination``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_contamination_CalculateContamination.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_contamination_CalculateContamination.php>`_
:Versions: 4.1.8.1, 4.1.7.0, 4.1.6.0, 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.7.0
:Authors: Hollizeck Sebastian
:Citations: TBD
:Created: 2019-09-09
:Updated: 2019-09-09


Outputs
-----------

=======  ======  =========================
name     type    documentation
=======  ======  =========================
contOut  File    contamination Table
segOut   File    segmentation based on baf
=======  ======  =========================


Additional configuration (inputs)
---------------------------------

====================  =======================  ==================================  ==========  =============================================================================================================================================
name                  type                     prefix                                position  documentation
====================  =======================  ==================================  ==========  =============================================================================================================================================
pileupTable           File                     -I                                              pileup table from summarize pileup
javaOptions           Optional<Array<String>>
compression_level     Optional<Integer>                                                        Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
contaminationTable    Optional<File>           --contamination-table                           Tables containing contamination information.
statsFile             Optional<File>           --stats                                         The Mutect stats file output by Mutect2
readOrientationModel  Optional<File>           --orientation-bias-artifact-priors              One or more .tar.gz files containing tables of prior artifact probabilities for the read orientation filter model, one table per tumor sample
segmentationFileOut   Optional<Filename>       --tumor-segmentation                            Reference sequence file
contaminationFileOut  Optional<Filename>       -O                                           2
====================  =======================  ==================================  ==========  =============================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4CalculateContamination {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       File? contaminationTable
       File? statsFile
       File? readOrientationModel
       File pileupTable
       String? segmentationFileOut
       String? contaminationFileOut
     }
     command <<<
       set -e
       gatk CalculateContamination \
         --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         ~{if defined(contaminationTable) then ("--contamination-table '" + contaminationTable + "'") else ""} \
         ~{if defined(statsFile) then ("--stats '" + statsFile + "'") else ""} \
         ~{if defined(readOrientationModel) then ("--orientation-bias-artifact-priors '" + readOrientationModel + "'") else ""} \
         -I '~{pileupTable}' \
         --tumor-segmentation '~{select_first([segmentationFileOut, "~{basename(pileupTable)}.mutect2_segments"])}' \
         -O '~{select_first([contaminationFileOut, "~{basename(pileupTable)}.mutect2_contamination"])}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.7.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File contOut = select_first([contaminationFileOut, "~{basename(pileupTable)}.mutect2_contamination"])
       File segOut = select_first([segmentationFileOut, "~{basename(pileupTable)}.mutect2_segments"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: CalculateContamination'
   doc: |-
     Calculates the fraction of reads coming from cross-sample contamination, given results from GetPileupSummaries. The resulting contamination table is used with FilterMutectCalls.

     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow. See Tutorial#11136 for a step-by-step description of the workflow and Article#11127 for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the Mutect2 WDL scripts directory.

     This tool borrows from ContEst by Cibulskis et al the idea of estimating contamination from ref reads at hom alt sites. However, ContEst uses a probabilistic model that assumes a diploid genotype with no copy number variation and independent contaminating reads. That is, ContEst assumes that each contaminating read is drawn randomly and independently from a different human. This tool uses a simpler estimate of contamination that relaxes these assumptions. In particular, it works in the presence of copy number variations and with an arbitrary number of contaminating samples. In addition, this tool is designed to work well with no matched normal data. However, one can run GetPileupSummaries on a matched normal bam file and input the result to this tool.

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
   - id: contaminationTable
     label: contaminationTable
     doc: Tables containing contamination information.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --contamination-table
   - id: statsFile
     label: statsFile
     doc: The Mutect stats file output by Mutect2
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --stats
   - id: readOrientationModel
     label: readOrientationModel
     doc: |-
       One or more .tar.gz files containing tables of prior artifact probabilities for the read orientation filter model, one table per tumor sample
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --orientation-bias-artifact-priors
   - id: pileupTable
     label: pileupTable
     doc: pileup table from summarize pileup
     type: File
     inputBinding:
       prefix: -I
   - id: segmentationFileOut
     label: segmentationFileOut
     doc: Reference sequence file
     type:
     - string
     - 'null'
     default: generated.mutect2_segments
     inputBinding:
       prefix: --tumor-segmentation
       valueFrom: $(inputs.pileupTable.basename).mutect2_segments
   - id: contaminationFileOut
     label: contaminationFileOut
     type:
     - string
     - 'null'
     default: generated.mutect2_contamination
     inputBinding:
       prefix: -O
       position: 2
       valueFrom: $(inputs.pileupTable.basename).mutect2_contamination

   outputs:
   - id: contOut
     label: contOut
     doc: contamination Table
     type: File
     outputBinding:
       glob: $(inputs.pileupTable.basename).mutect2_contamination
       loadContents: false
   - id: segOut
     label: segOut
     doc: segmentation based on baf
     type: File
     outputBinding:
       glob: $(inputs.pileupTable.basename).mutect2_segments
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - CalculateContamination
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4CalculateContamination


