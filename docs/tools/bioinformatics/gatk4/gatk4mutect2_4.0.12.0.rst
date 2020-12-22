:orphan:

GATK4: MuTect2
=============================

``Gatk4Mutect2`` · *1 contributor · 7 versions*

Call somatic short variants via local assembly of haplotypes. Short variants include single nucleotide (SNV)
and insertion and deletion (indel) variants. The caller combines the DREAM challenge-winning somatic
genotyping engine of the original MuTect (Cibulskis et al., 2013) with the assembly-based machinery of HaplotypeCaller.

This tool is featured in the Somatic Short Mutation calling Best Practice Workflow. See Tutorial#11136
for a step-by-step description of the workflow and Article#11127 for an overview of what traditional
somatic calling entails. For the latest pipeline scripts, see the Mutect2 WDL scripts directory.
Although we present the tool for somatic calling, it may apply to other contexts,
such as mitochondrial variant calling.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.mutect2.versions import GatkMutect2_4_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4mutect2_step",
           GatkMutect2_4_0(
               tumor=None,
               tumorName=None,
               normal=None,
               normalName=None,
               reference=None,
           )
       )
       wf.output("out", source=gatk4mutect2_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4Mutect2:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4Mutect2 > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       normal: normal.bam
       normalName: <value>
       reference: reference.fasta
       tumor: tumor.bam
       tumorName: <value>




5. Run Gatk4Mutect2 with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4Mutect2





Information
------------

:ID: ``Gatk4Mutect2``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.10.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.10.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php>`_
:Versions: 4.1.8.1, 4.1.7.0, 4.1.6.0, 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.0.12.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24


Outputs
-----------

======  ============  =================
name    type          documentation
======  ============  =================
out     Gzipped<VCF>  To determine type
======  ============  =================


Additional configuration (inputs)
---------------------------------

========================  =======================  ===============================  ==========  ==============================================================================================================================================================
name                      type                     prefix                             position  documentation
========================  =======================  ===============================  ==========  ==============================================================================================================================================================
tumor                     IndexedBam               -I                                        6  BAM/SAM/CRAM file containing reads
tumorName                 String                   -tumor                                    6  BAM sample name of tumor. May be URL-encoded as output by GetSampleName with -encode.
normal                    IndexedBam               -I                                        5  BAM/SAM/CRAM file containing reads
normalName                String                   -normal                                   6  BAM sample name of normal. May be URL-encoded as output by GetSampleName with -encode.
reference                 FastaWithIndexes         -R                                        8  Reference sequence file
javaOptions               Optional<Array<String>>
compression_level         Optional<Integer>                                                     Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
intervals                 Optional<bed>            -L                                        7  One or more genomic intervals over which to operate
outputFilename            Optional<Filename>       -O                                       20
germlineResource          Optional<IndexedVCF>     --germline-resource                      10
afOfAllelesNotInResource  Optional<Float>          --af-of-alleles-not-in-resource          11  Population allele fraction assigned to alleles not found in germline resource. Please see docs/mutect/mutect2.pdf fora derivation of the default value.
panelOfNormals            Optional<IndexedVCF>     --panel-of-normals                       10  A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants.
========================  =======================  ===============================  ==========  ==============================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4Mutect2 {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       File tumor
       File tumor_bai
       String tumorName
       File normal
       File normal_bai
       String normalName
       File? intervals
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       String? outputFilename
       File? germlineResource
       File? germlineResource_idx
       Float? afOfAllelesNotInResource
       File? panelOfNormals
       File? panelOfNormals_idx
     }
     command <<<
       set -e
       gatk Mutect2 \
         --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         -I '~{normal}' \
         -I '~{tumor}' \
         -tumor '~{tumorName}' \
         -normal '~{normalName}' \
         ~{if defined(intervals) then ("-L '" + intervals + "'") else ""} \
         -R '~{reference}' \
         ~{if defined(germlineResource) then ("--germline-resource '" + germlineResource + "'") else ""} \
         ~{if defined(panelOfNormals) then ("--panel-of-normals '" + panelOfNormals + "'") else ""} \
         ~{if defined(afOfAllelesNotInResource) then ("--af-of-alleles-not-in-resource " + afOfAllelesNotInResource) else ''} \
         -O '~{select_first([outputFilename, "generated.vcf.gz"])}'
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
       File out = select_first([outputFilename, "generated.vcf.gz"])
       File out_tbi = select_first([outputFilename, "generated.vcf.gz"]) + ".tbi"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: MuTect2'
   doc: |-
     Call somatic short variants via local assembly of haplotypes. Short variants include single nucleotide (SNV)
     and insertion and deletion (indel) variants. The caller combines the DREAM challenge-winning somatic
     genotyping engine of the original MuTect (Cibulskis et al., 2013) with the assembly-based machinery of HaplotypeCaller.

     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow. See Tutorial#11136
     for a step-by-step description of the workflow and Article#11127 for an overview of what traditional
     somatic calling entails. For the latest pipeline scripts, see the Mutect2 WDL scripts directory.
     Although we present the tool for somatic calling, it may apply to other contexts,
     such as mitochondrial variant calling.

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
   - id: tumor
     label: tumor
     doc: BAM/SAM/CRAM file containing reads
     type: File
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       prefix: -I
       position: 6
   - id: tumorName
     label: tumorName
     doc: |-
       BAM sample name of tumor. May be URL-encoded as output by GetSampleName with -encode.
     type: string
     inputBinding:
       prefix: -tumor
       position: 6
   - id: normal
     label: normal
     doc: BAM/SAM/CRAM file containing reads
     type: File
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       prefix: -I
       position: 5
   - id: normalName
     label: normalName
     doc: |-
       BAM sample name of normal. May be URL-encoded as output by GetSampleName with -encode.
     type: string
     inputBinding:
       prefix: -normal
       position: 6
   - id: intervals
     label: intervals
     doc: One or more genomic intervals over which to operate
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -L
       position: 7
   - id: reference
     label: reference
     doc: Reference sequence file
     type: File
     secondaryFiles:
     - pattern: .fai
     - pattern: .amb
     - pattern: .ann
     - pattern: .bwt
     - pattern: .pac
     - pattern: .sa
     - pattern: ^.dict
     inputBinding:
       prefix: -R
       position: 8
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.vcf.gz
     inputBinding:
       prefix: -O
       position: 20
   - id: germlineResource
     label: germlineResource
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .idx
     inputBinding:
       prefix: --germline-resource
       position: 10
   - id: afOfAllelesNotInResource
     label: afOfAllelesNotInResource
     doc: |-
       Population allele fraction assigned to alleles not found in germline resource. Please see docs/mutect/mutect2.pdf fora derivation of the default value.
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --af-of-alleles-not-in-resource
       position: 11
   - id: panelOfNormals
     label: panelOfNormals
     doc: |-
       A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants.
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .idx
     inputBinding:
       prefix: --panel-of-normals
       position: 10

   outputs:
   - id: out
     label: out
     doc: To determine type
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: generated.vcf.gz
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - Mutect2
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4Mutect2


