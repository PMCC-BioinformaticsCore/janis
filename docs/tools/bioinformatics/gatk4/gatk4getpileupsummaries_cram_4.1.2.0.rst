:orphan:

GATK4: GetPileupSummaries
========================================================

``Gatk4GetPileupSummaries_cram`` · *1 contributor · 6 versions*

Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results can be used with CalculateContamination.
The tool requires a common germline variant sites VCF, e.g. the gnomAD resource, with population allele frequencies (AF) in the INFO field. This resource must contain only biallelic SNPs and can be an eight-column sites-only VCF. The tool ignores the filter status of the sites. See the GATK Resource Bundle for an example human file.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.getpileupsummaries.versions import Gatk4GetPileUpSummariesCram_4_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4getpileupsummaries_cram_step",
           Gatk4GetPileUpSummariesCram_4_1_2(
               bam=None,
               sites=None,
           )
       )
       wf.output("out", source=gatk4getpileupsummaries_cram_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4GetPileupSummaries_cram:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4GetPileupSummaries_cram > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam:
       - bam_0.cram
       - bam_1.cram
       sites: sites.vcf.gz




5. Run Gatk4GetPileupSummaries_cram with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4GetPileupSummaries_cram





Information
------------

:ID: ``Gatk4GetPileupSummaries_cram``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/org_broadinstitute_hellbender_tools_walkers_contamination_GetPileupSummaries.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/org_broadinstitute_hellbender_tools_walkers_contamination_GetPileupSummaries.php>`_
:Versions: 4.1.8.1, 4.1.7.0, 4.1.6.0, 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.2.0
:Authors: Hollizeck Sebastian
:Citations: TBD
:Created: 2019-09-09
:Updated: 2019-09-09


Outputs
-----------

======  ========  ================================
name    type      documentation
======  ========  ================================
out     TextFile  Table containing the pileup info
======  ========  ================================


Additional configuration (inputs)
---------------------------------

=================  ==========================  ===========  ==========  ========================================================================================
name               type                        prefix         position  documentation
=================  ==========================  ===========  ==========  ========================================================================================
bam                Array<CramPair>             -I                    0  The SAM/BAM/CRAM file containing reads.
sites              Gzipped<VCF>                -V                       sites of common biallelic variants
javaOptions        Optional<Array<String>>
compression_level  Optional<Integer>                                    Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
sampleName         Optional<String>                                     Used for naming purposes
intervals          Optional<bed>               --intervals              -L (BASE) One or more genomic intervals over which to operate
pileupTableOut     Optional<Filename>          -O                    1
reference          Optional<FastaWithIndexes>  -R                       reference to use when decoding CRAMS
=================  ==========================  ===========  ==========  ========================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4GetPileupSummaries_cram {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       Array[File] bam
       Array[File] bam_crai
       String? sampleName
       File sites
       File sites_tbi
       File? intervals
       String? pileupTableOut
       File? reference
       File? reference_fai
       File? reference_amb
       File? reference_ann
       File? reference_bwt
       File? reference_pac
       File? reference_sa
       File? reference_dict
     }
     command <<<
       set -e
       gatk GetPileupSummaries \
         --java-options '-Xmx~{((select_first([runtime_memory, 64, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         ~{if length(bam) > 0 then "-I '" + sep("' -I '", bam) + "'" else ""} \
         -V '~{sites}' \
         ~{if defined(intervals) then ("--intervals '" + intervals + "'") else ""} \
         ~{if defined(reference) then ("-R '" + reference + "'") else ""} \
         -O '~{select_first([pileupTableOut, "~{sep(".", select_all([select_first([sampleName, "generated"])]))}.txt"])}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.2.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 64, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([pileupTableOut, "~{sep(".", select_all([select_first([sampleName, "generated"])]))}.txt"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: GetPileupSummaries'
   doc: |-
     Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results can be used with CalculateContamination.
     The tool requires a common germline variant sites VCF, e.g. the gnomAD resource, with population allele frequencies (AF) in the INFO field. This resource must contain only biallelic SNPs and can be an eight-column sites-only VCF. The tool ignores the filter status of the sites. See the GATK Resource Bundle for an example human file.

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
   - id: bam
     label: bam
     doc: The SAM/BAM/CRAM file containing reads.
     type:
       type: array
       inputBinding:
         prefix: -I
       items: File
     inputBinding:
       position: 0
   - id: sampleName
     label: sampleName
     doc: Used for naming purposes
     type:
     - string
     - 'null'
   - id: sites
     label: sites
     doc: sites of common biallelic variants
     type: File
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       prefix: -V
   - id: intervals
     label: intervals
     doc: -L (BASE) One or more genomic intervals over which to operate
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --intervals
   - id: pileupTableOut
     label: pileupTableOut
     type:
     - string
     - 'null'
     default: generated.txt
     inputBinding:
       prefix: -O
       position: 1
       valueFrom: |-
         $([[inputs.sampleName, "generated"].filter(function (inner) { return inner != null })[0]].filter(function (inner) { return inner != null }).join(".")).txt
   - id: reference
     label: reference
     doc: reference to use when decoding CRAMS
     type:
     - File
     - 'null'
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

   outputs:
   - id: out
     label: out
     doc: Table containing the pileup info
     type: File
     outputBinding:
       glob: |-
         $([[inputs.sampleName, "generated"].filter(function (inner) { return inner != null })[0]].filter(function (inner) { return inner != null }).join(".")).txt
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - GetPileupSummaries
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 64, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4GetPileupSummaries_cram


