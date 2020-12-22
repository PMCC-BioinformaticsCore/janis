:orphan:

Manta
==================

``manta_cram`` · *1 contributor · 2 versions*

Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads.
It is optimized for analysis of germline variation in small sets of individuals and somatic
variation in tumor/normal sample pairs. Manta discovers, assembles and scores large-scale SVs,
medium-sized indels and large insertions within a single efficient workflow. The method is
designed for rapid analysis on standard compute hardware: NA12878 at 50x genomic coverage is
analyzed in less than 20 minutes on a 20 core server, and most WGS tumor/normal analyses
can be completed within 2 hours. Manta combines paired and split-read evidence during SV
discovery and scoring to improve accuracy, but does not require split-reads or successful
breakpoint assemblies to report a variant in cases where there is strong evidence otherwise.

It provides scoring models for germline variants in small sets of diploid samples and somatic
variants in matched tumor/normal sample pairs. There is experimental support for analysis of
unmatched tumor samples as well. Manta accepts input read mappings from BAM or CRAM files and
reports all SV and indel inferences in VCF 4.1 format. See the user guide for a full description
of capabilities and limitations.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.illumina.manta.manta import MantaCram_1_5_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "manta_cram_step",
           MantaCram_1_5_0(
               bam=None,
               reference=None,
           )
       )
       wf.output("python", source=manta_cram_step.python)
       wf.output("pickle", source=manta_cram_step.pickle)
       wf.output("candidateSV", source=manta_cram_step.candidateSV)
       wf.output("candidateSmallIndels", source=manta_cram_step.candidateSmallIndels)
       wf.output("diploidSV", source=manta_cram_step.diploidSV)
       wf.output("alignmentStatsSummary", source=manta_cram_step.alignmentStatsSummary)
       wf.output("svCandidateGenerationStats", source=manta_cram_step.svCandidateGenerationStats)
       wf.output("svLocusGraphStats", source=manta_cram_step.svLocusGraphStats)
       wf.output("somaticSVs", source=manta_cram_step.somaticSVs)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for manta_cram:

.. code-block:: bash

   # user inputs
   janis inputs manta_cram > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.cram
       reference: reference.fasta




5. Run manta_cram with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       manta_cram





Information
------------

:ID: ``manta_cram``
:URL: `https://github.com/Illumina/manta <https://github.com/Illumina/manta>`_
:Versions: 1.5.0, 1.4.0
:Container: michaelfranklin/manta:1.5.0
:Authors: Michael Franklin
:Citations: Chen, X. et al. (2016) Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics, 32, 1220-1222. doi:10.1093/bioinformatics/btv710
:DOI:  doi:10.1093/bioinformatics/btv710
:Created: 2019-02-12
:Updated: 2019-02-19


Outputs
-----------

==========================  ======================  ===============
name                        type                    documentation
==========================  ======================  ===============
python                      File
pickle                      File
candidateSV                 Gzipped<VCF>
candidateSmallIndels        Gzipped<VCF>
diploidSV                   Gzipped<VCF>
alignmentStatsSummary       File
svCandidateGenerationStats  tsv
svLocusGraphStats           tsv
somaticSVs                  Optional<Gzipped<VCF>>
==========================  ======================  ===============


Additional configuration (inputs)
---------------------------------

==============  ======================  ================  ==========  ====================================================================================================================================================================================================================================================================================================================================================
name            type                    prefix              position  documentation
==============  ======================  ================  ==========  ====================================================================================================================================================================================================================================================================================================================================================
bam             CramPair                --bam                      1  FILE Normal sample BAM or CRAM file. May be specified more than once, multiple inputs will be treated as each BAM file representing a different sample. [optional] (no default)
reference       FastaFai                --referenceFasta           1  samtools-indexed reference fasta file [required]
config          Optional<File>          --config                   1  provide a configuration file to override defaults in global config file (/opt/conda/share/manta-1.2.1-0/bin/configManta.py.ini)
runDir          Optional<Filename>      --runDir                   1  Run script and run output will be written to this directory [required] (default: MantaWorkflow)
tumorBam        Optional<CramPair>      --tumorBam                 1  Tumor sample BAM or CRAM file. Only up to one tumor bam file accepted. [optional=null]
exome           Optional<Boolean>       --exome                    1  Set options for WES input: turn off depth filters
rna             Optional<Boolean>       --rna                      1  Set options for RNA-Seq input. Must specify exactly one bam input file
unstrandedRNA   Optional<Boolean>       --unstrandedRNA            1  Set if RNA-Seq input is unstranded: Allows splice-junctions on either strand
outputContig    Optional<Boolean>       --outputContig             1  Output assembled contig sequences in VCF file
callRegions     Optional<Gzipped<bed>>  --callRegions              1  Optionally provide a bgzip-compressed/tabix-indexed BED file containing the set of regions to call. No VCF output will be provided outside of these regions. The full genome will still be used to estimate statistics from the input (such as expected depth per chromosome). Only one BED file may be specified. (default: call the entire genome)
mode            Optional<String>        --mode                     3  (-m) select run mode (local|sge)
quiet           Optional<Boolean>       --quiet                    3  Don't write any log output to stderr (but still write to workspace/pyflow.data/logs/pyflow_log.txt)
queue           Optional<String>        --queue                    3  (-q) specify scheduler queue name
memgb           Optional<Integer>       --memGb                    3  (-g) gigabytes of memory available to run workflow -- only meaningful in local mode, must be an integer (default: Estimate the total memory for this node for local  mode, 'unlimited' for sge mode)
maxTaskRuntime  Optional<String>        --maxTaskRuntime           3  (format: hh:mm:ss) Specify scheduler max runtime per task, argument is provided to the 'h_rt' resource limit if using SGE (no default)
==============  ======================  ================  ==========  ====================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task manta_cram {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File? config
       File bam
       File bam_crai
       String? runDir
       File reference
       File reference_fai
       File? tumorBam
       File? tumorBam_crai
       Boolean? exome
       Boolean? rna
       Boolean? unstrandedRNA
       Boolean? outputContig
       File? callRegions
       File? callRegions_tbi
       String? mode
       Boolean? quiet
       String? queue
       Int? memgb
       String? maxTaskRuntime
     }
     command <<<
       set -e
        \
         configManta.py \
         ~{if defined(config) then ("--config " + config) else ''} \
         --bam ~{bam} \
         --runDir ~{select_first([runDir, "generated"])} \
         --referenceFasta ~{reference} \
         ~{if defined(tumorBam) then ("--tumorBam " + tumorBam) else ''} \
         ~{if (defined(exome) && select_first([exome])) then "--exome" else ""} \
         ~{if (defined(rna) && select_first([rna])) then "--rna" else ""} \
         ~{if (defined(unstrandedRNA) && select_first([unstrandedRNA])) then "--unstrandedRNA" else ""} \
         ~{if (defined(outputContig) && select_first([outputContig])) then "--outputContig" else ""} \
         ~{if defined(callRegions) then ("--callRegions " + callRegions) else ''} \
         ;~{select_first([runDir, "generated"])}/runWorkflow.py \
         ~{if defined(select_first([mode, "local"])) then ("--mode " + select_first([mode, "local"])) else ''} \
         ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
         ~{if defined(queue) then ("--queue " + queue) else ''} \
         ~{if defined(memgb) then ("--memGb " + memgb) else ''} \
         ~{if defined(maxTaskRuntime) then ("--maxTaskRuntime " + maxTaskRuntime) else ''} \
         -j ~{select_first([runtime_cpu, 4])}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "michaelfranklin/manta:1.5.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4, 4])}G"
       preemptible: 2
     }
     output {
       File python = (select_first([runDir, "generated"]) + "/runWorkflow.py")
       File pickle = (select_first([runDir, "generated"]) + "/runWorkflow.py.config.pickle")
       File candidateSV = (select_first([runDir, "generated"]) + "/results/variants/candidateSV.vcf.gz")
       File candidateSV_tbi = (select_first([runDir, "generated"]) + "/results/variants/candidateSV.vcf.gz") + ".tbi"
       File candidateSmallIndels = (select_first([runDir, "generated"]) + "/results/variants/candidateSmallIndels.vcf.gz")
       File candidateSmallIndels_tbi = (select_first([runDir, "generated"]) + "/results/variants/candidateSmallIndels.vcf.gz") + ".tbi"
       File diploidSV = (select_first([runDir, "generated"]) + "/results/variants/diploidSV.vcf.gz")
       File diploidSV_tbi = (select_first([runDir, "generated"]) + "/results/variants/diploidSV.vcf.gz") + ".tbi"
       File alignmentStatsSummary = (select_first([runDir, "generated"]) + "/results/stats/alignmentStatsSummary.txt")
       File svCandidateGenerationStats = (select_first([runDir, "generated"]) + "/results/stats/svCandidateGenerationStats.tsv")
       File svLocusGraphStats = (select_first([runDir, "generated"]) + "/results/stats/svLocusGraphStats.tsv")
       File? somaticSVs = (select_first([runDir, "generated"]) + "/results/variants/somaticSV.vcf.gz")
       File? somaticSVs_tbi = if defined((select_first([runDir, "generated"]) + "/results/variants/somaticSV.vcf.gz")) then ((select_first([runDir, "generated"]) + "/results/variants/somaticSV.vcf.gz") + ".tbi") else None
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Manta
   doc: |-
     Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads.
     It is optimized for analysis of germline variation in small sets of individuals and somatic
     variation in tumor/normal sample pairs. Manta discovers, assembles and scores large-scale SVs,
     medium-sized indels and large insertions within a single efficient workflow. The method is
     designed for rapid analysis on standard compute hardware: NA12878 at 50x genomic coverage is
     analyzed in less than 20 minutes on a 20 core server, and most WGS tumor/normal analyses
     can be completed within 2 hours. Manta combines paired and split-read evidence during SV
     discovery and scoring to improve accuracy, but does not require split-reads or successful
     breakpoint assemblies to report a variant in cases where there is strong evidence otherwise.

     It provides scoring models for germline variants in small sets of diploid samples and somatic
     variants in matched tumor/normal sample pairs. There is experimental support for analysis of
     unmatched tumor samples as well. Manta accepts input read mappings from BAM or CRAM files and
     reports all SV and indel inferences in VCF 4.1 format. See the user guide for a full description
     of capabilities and limitations.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: michaelfranklin/manta:1.5.0

   inputs:
   - id: config
     label: config
     doc: |-
       provide a configuration file to override defaults in global config file (/opt/conda/share/manta-1.2.1-0/bin/configManta.py.ini)
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --config
       position: 1
       shellQuote: false
   - id: bam
     label: bam
     doc: |-
       FILE Normal sample BAM or CRAM file. May be specified more than once, multiple inputs will be treated as each BAM file representing a different sample. [optional] (no default)
     type: File
     secondaryFiles:
     - pattern: .crai
     inputBinding:
       prefix: --bam
       position: 1
       shellQuote: false
   - id: runDir
     label: runDir
     doc: |-
       Run script and run output will be written to this directory [required] (default: MantaWorkflow)
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --runDir
       position: 1
       shellQuote: false
   - id: reference
     label: reference
     doc: samtools-indexed reference fasta file [required]
     type: File
     secondaryFiles:
     - pattern: .fai
     inputBinding:
       prefix: --referenceFasta
       position: 1
       shellQuote: false
   - id: tumorBam
     label: tumorBam
     doc: |-
       Tumor sample BAM or CRAM file. Only up to one tumor bam file accepted. [optional=null]
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .crai
     inputBinding:
       prefix: --tumorBam
       position: 1
       shellQuote: false
   - id: exome
     label: exome
     doc: 'Set options for WES input: turn off depth filters'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --exome
       position: 1
       shellQuote: false
   - id: rna
     label: rna
     doc: Set options for RNA-Seq input. Must specify exactly one bam input file
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --rna
       position: 1
       shellQuote: false
   - id: unstrandedRNA
     label: unstrandedRNA
     doc: 'Set if RNA-Seq input is unstranded: Allows splice-junctions on either strand'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --unstrandedRNA
       position: 1
       shellQuote: false
   - id: outputContig
     label: outputContig
     doc: Output assembled contig sequences in VCF file
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --outputContig
       position: 1
       shellQuote: false
   - id: callRegions
     label: callRegions
     doc: |-
       Optionally provide a bgzip-compressed/tabix-indexed BED file containing the set of regions to call. No VCF output will be provided outside of these regions. The full genome will still be used to estimate statistics from the input (such as expected depth per chromosome). Only one BED file may be specified. (default: call the entire genome)
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       prefix: --callRegions
       position: 1
       shellQuote: false
   - id: mode
     label: mode
     doc: (-m) select run mode (local|sge)
     type: string
     default: local
     inputBinding:
       prefix: --mode
       position: 3
       shellQuote: false
   - id: quiet
     label: quiet
     doc: |-
       Don't write any log output to stderr (but still write to workspace/pyflow.data/logs/pyflow_log.txt)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --quiet
       position: 3
       shellQuote: false
   - id: queue
     label: queue
     doc: (-q) specify scheduler queue name
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --queue
       position: 3
       shellQuote: false
   - id: memgb
     label: memgb
     doc: |-
       (-g) gigabytes of memory available to run workflow -- only meaningful in local mode, must be an integer (default: Estimate the total memory for this node for local  mode, 'unlimited' for sge mode)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --memGb
       position: 3
       shellQuote: false
   - id: maxTaskRuntime
     label: maxTaskRuntime
     doc: |-
       (format: hh:mm:ss) Specify scheduler max runtime per task, argument is provided to the 'h_rt' resource limit if using SGE (no default)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --maxTaskRuntime
       position: 3
       shellQuote: false

   outputs:
   - id: python
     label: python
     type: File
     outputBinding:
       glob: $((inputs.runDir + "/runWorkflow.py"))
       outputEval: $((inputs.runDir.basename + "/runWorkflow.py"))
       loadContents: false
   - id: pickle
     label: pickle
     type: File
     outputBinding:
       glob: $((inputs.runDir + "/runWorkflow.py.config.pickle"))
       outputEval: $((inputs.runDir.basename + "/runWorkflow.py.config.pickle"))
       loadContents: false
   - id: candidateSV
     label: candidateSV
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $((inputs.runDir + "/results/variants/candidateSV.vcf.gz"))
       outputEval: $((inputs.runDir.basename + "/results/variants/candidateSV.vcf.gz"))
       loadContents: false
   - id: candidateSmallIndels
     label: candidateSmallIndels
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $((inputs.runDir + "/results/variants/candidateSmallIndels.vcf.gz"))
       outputEval: $((inputs.runDir.basename + "/results/variants/candidateSmallIndels.vcf.gz"))
       loadContents: false
   - id: diploidSV
     label: diploidSV
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $((inputs.runDir + "/results/variants/diploidSV.vcf.gz"))
       outputEval: $((inputs.runDir.basename + "/results/variants/diploidSV.vcf.gz"))
       loadContents: false
   - id: alignmentStatsSummary
     label: alignmentStatsSummary
     type: File
     outputBinding:
       glob: $((inputs.runDir + "/results/stats/alignmentStatsSummary.txt"))
       outputEval: $((inputs.runDir.basename + "/results/stats/alignmentStatsSummary.txt"))
       loadContents: false
   - id: svCandidateGenerationStats
     label: svCandidateGenerationStats
     type: File
     outputBinding:
       glob: $((inputs.runDir + "/results/stats/svCandidateGenerationStats.tsv"))
       outputEval: $((inputs.runDir.basename + "/results/stats/svCandidateGenerationStats.tsv"))
       loadContents: false
   - id: svLocusGraphStats
     label: svLocusGraphStats
     type: File
     outputBinding:
       glob: $((inputs.runDir + "/results/stats/svLocusGraphStats.tsv"))
       outputEval: $((inputs.runDir.basename + "/results/stats/svLocusGraphStats.tsv"))
       loadContents: false
   - id: somaticSVs
     label: somaticSVs
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $((inputs.runDir + "/results/variants/somaticSV.vcf.gz"))
       outputEval: $((inputs.runDir.basename + "/results/variants/somaticSV.vcf.gz"))
       loadContents: false
   stdout: _stdout
   stderr: _stderr
   arguments:
   - position: 0
     valueFrom: configManta.py
     shellQuote: false
   - position: 2
     valueFrom: $(";{runDir}/runWorkflow.py".replace(/\{runDir\}/g, inputs.runDir))
     shellQuote: false
   - prefix: -j
     position: 3
     valueFrom: $([inputs.runtime_cpu, 4].filter(function (inner) { return inner != null
       })[0])
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: manta_cram


