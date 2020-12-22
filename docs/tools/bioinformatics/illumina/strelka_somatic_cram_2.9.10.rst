:orphan:

Strelka (Somatic)
========================================

``strelka_somatic_cram`` · *1 contributor · 2 versions*

Usage: configureStrelkaSomaticWorkflow.py [options]
Version: 2.9.10
This script configures Strelka somatic small variant calling.
You must specify an alignment file (BAM or CRAM) for each sample of a matched tumor-normal pair.
Configuration will produce a workflow run script which can execute the workflow on a single node or through
sge and resume any interrupted execution.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.illumina.strelkasomatic.strelkasomatic import StrelkaSomaticCram_2_9_10

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "strelka_somatic_cram_step",
           StrelkaSomaticCram_2_9_10(
               normalBam=None,
               tumorBam=None,
               reference=None,
           )
       )
       wf.output("configPickle", source=strelka_somatic_cram_step.configPickle)
       wf.output("script", source=strelka_somatic_cram_step.script)
       wf.output("stats", source=strelka_somatic_cram_step.stats)
       wf.output("indels", source=strelka_somatic_cram_step.indels)
       wf.output("snvs", source=strelka_somatic_cram_step.snvs)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for strelka_somatic_cram:

.. code-block:: bash

   # user inputs
   janis inputs strelka_somatic_cram > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       normalBam: normalBam.cram
       reference: reference.fasta
       tumorBam: tumorBam.cram




5. Run strelka_somatic_cram with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       strelka_somatic_cram





Information
------------

:ID: ``strelka_somatic_cram``
:URL: *No URL to the documentation was provided*
:Versions: 2.9.10, 2.9.9
:Container: michaelfranklin/strelka:2.9.10
:Authors: Michael Franklin
:Citations: None
:Created: 2019-05-27
:Updated: 2019-10-14


Outputs
-----------

============  ============  ===========================================================================================================================================================================================================================================
name          type          documentation
============  ============  ===========================================================================================================================================================================================================================================
configPickle  File
script        File
stats         tsv           A tab-delimited report of various internal statistics from the variant calling process: Runtime information accumulated for each genome segment, excluding auxiliary steps such as BAM indexing and vcf merging. Indel candidacy statistics
indels        Gzipped<VCF>
snvs          Gzipped<VCF>
============  ============  ===========================================================================================================================================================================================================================================


Additional configuration (inputs)
---------------------------------

=====================  =============================  ========================  ==========  ====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                   type                           prefix                      position  documentation
=====================  =============================  ========================  ==========  ====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
normalBam              CramPair                       --normalBam=                       1  Normal sample BAM or CRAM file. (no default)
tumorBam               CramPair                       --tumourBam=                       1  (--tumorBam)  Tumor sample BAM or CRAM file. [required] (no default)
reference              FastaFai                       --referenceFasta=                  1  samtools-indexed reference fasta file [required]
rundir                 Optional<Filename>             --runDir=                          1  Name of directory to be created where all workflow scripts and output will be written. Each analysis requires a separate directory. (default: StrelkaSomaticWorkflow)
region                 Optional<Array<String>>        --region                           1  Limit the analysis to one or more genome region(s) for debugging purposes. If this argument is provided multiple times the union of all specified regions will be analyzed. All regions must be non-overlapping to get a meaningful result. Examples: '--region chr20' (whole chromosome), '--region chr2:100-2000 --region chr3:2500-3000' (two regions)'. If this option is specified (one or more times) together with the 'callRegions' BED file,then all region arguments will be intersected with the callRegions BED track.
config                 Optional<File>                 --config=                          1  provide a configuration file to override defaults in global config file (/opt/strelka/bin/configureStrelkaSomaticWorkflow.py.ini)
outputcallableregions  Optional<Boolean>              --outputCallableRegions            1  Output a bed file describing somatic callable regions of the genome
indelCandidates        Optional<Array<Gzipped<VCF>>>  --indelCandidates=                 1  Specify a VCF of candidate indel alleles. These alleles are always evaluated but only reported in the output when they are inferred to exist in the sample. The VCF must be tabix indexed. All indel alleles must be left-shifted/normalized, any unnormalized alleles will be ignored. This option may be specified more than once, multiple input VCFs will be merged. (default: None)
forcedgt               Optional<Array<Gzipped<VCF>>>  --forcedGT=                        1  Specify a VCF of candidate alleles. These alleles are always evaluated and reported even if they are unlikely to exist in the sample. The VCF must be tabix indexed. All indel alleles must be left- shifted/normalized, any unnormalized allele will trigger a runtime error. This option may be specified more than once, multiple input VCFs will be merged. Note that for any SNVs provided in the VCF, the SNV site will be reported (and for gVCF, excluded from block compression), but the specific SNV alleles are ignored. (default: None)
targeted               Optional<Boolean>              --targeted                         1  Set options for other targeted input: note in particular that this flag turns off high-depth filters
exome                  Optional<Boolean>              --exome                            1  Set options for exome: note in particular that this flag turns off high-depth filters
callRegions            Optional<Gzipped<bed>>         --callRegions=                     1  Optionally provide a bgzip-compressed/tabix-indexed BED file containing the set of regions to call. No VCF output will be provided outside of these regions. The full genome will still be used to estimate statistics from the input (such as expected depth per chromosome). Only one BED file may be specified. (default: call the entire genome)
noisevcf               Optional<Gzipped<VCF>>         --noiseVcf=                        1  Noise vcf file (submit argument multiple times for more than one file)
scansizemb             Optional<Integer>              --scanSizeMb=                      1  Maximum sequence region size (in megabases) scanned by each task during genome variant calling. (default: 12)
callmemmb              Optional<Integer>              --callMemMb=                       1  Set variant calling task memory limit (in megabytes). It is not recommended to change the default in most cases, but this might be required for a sample of unusual depth.
retaintempfiles        Optional<Boolean>              --retainTempFiles                  1  Keep all temporary files (for workflow debugging)
disableevs             Optional<Boolean>              --disableEVS                       1  Disable empirical variant scoring (EVS).
reportevsfeatures      Optional<Boolean>              --reportEVSFeatures                1  Report all empirical variant scoring features in VCF output.
snvscoringmodelfile    Optional<File>                 --snvScoringModelFile=             1  Provide a custom empirical scoring model file for SNVs (default: /opt/strelka/share/config/somaticSNVScoringM odels.json)
indelscoringmodelfile  Optional<File>                 --indelScoringModelFile=           1  Provide a custom empirical scoring model file for indels (default: /opt/strelka/share/config/somaticInde lScoringModels.json)
mode                   Optional<String>               --mode                             3  (-m MODE)  select run mode (local|sge)
queue                  Optional<String>               --queue                            3  (-q QUEUE) specify scheduler queue name
memGb                  Optional<String>               --memGb                            3  (-g MEMGB) gigabytes of memory available to run workflow -- only meaningful in local mode, must be an integer (default: Estimate the total memory for this node for local mode, 'unlimited' for sge mode)
quiet                  Optional<Boolean>              --quiet                            3  Don't write any log output to stderr (but still write to workspace/pyflow.data/logs/pyflow_log.txt)
=====================  =============================  ========================  ==========  ====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task strelka_somatic_cram {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File normalBam
       File normalBam_crai
       File tumorBam
       File tumorBam_crai
       File reference
       File reference_fai
       String? rundir
       Array[String]? region
       File? config
       Boolean? outputcallableregions
       Array[File]? indelCandidates
       Array[File]? indelCandidates_tbi
       Array[File]? forcedgt
       Array[File]? forcedgt_tbi
       Boolean? targeted
       Boolean? exome
       File? callRegions
       File? callRegions_tbi
       File? noisevcf
       File? noisevcf_tbi
       Int? scansizemb
       Int? callmemmb
       Boolean? retaintempfiles
       Boolean? disableevs
       Boolean? reportevsfeatures
       File? snvscoringmodelfile
       File? indelscoringmodelfile
       String? mode
       String? queue
       String? memGb
       Boolean? quiet
     }
     command <<<
       set -e
        \
         'configureStrelkaSomaticWorkflow.py' \
         --normalBam='~{normalBam}' \
         --tumourBam='~{tumorBam}' \
         --referenceFasta='~{reference}' \
         --runDir='~{select_first([rundir, "generated"])}' \
         ~{if (defined(region) && length(select_first([region])) > 0) then "--region '" + sep("' --region '", select_first([region])) + "'" else ""} \
         ~{if defined(config) then ("--config='" + config + "'") else ""} \
         ~{if (defined(outputcallableregions) && select_first([outputcallableregions])) then "--outputCallableRegions" else ""} \
         ~{if (defined(indelCandidates) && length(select_first([indelCandidates])) > 0) then "--indelCandidates='" + sep("' --indelCandidates='", select_first([indelCandidates])) + "'" else ""} \
         ~{if (defined(forcedgt) && length(select_first([forcedgt])) > 0) then "--forcedGT='" + sep("' --forcedGT='", select_first([forcedgt])) + "'" else ""} \
         ~{if (defined(targeted) && select_first([targeted])) then "--targeted" else ""} \
         ~{if (defined(exome) && select_first([exome])) then "--exome" else ""} \
         ~{if defined(callRegions) then ("--callRegions='" + callRegions + "'") else ""} \
         ~{if defined(noisevcf) then ("--noiseVcf='" + noisevcf + "'") else ""} \
         ~{if defined(scansizemb) then ("--scanSizeMb=" + scansizemb) else ''} \
         ~{if defined(callmemmb) then ("--callMemMb=" + callmemmb) else ''} \
         ~{if select_first([retaintempfiles, false]) then "--retainTempFiles" else ""} \
         ~{if (defined(disableevs) && select_first([disableevs])) then "--disableEVS" else ""} \
         ~{if (defined(reportevsfeatures) && select_first([reportevsfeatures])) then "--reportEVSFeatures" else ""} \
         ~{if defined(snvscoringmodelfile) then ("--snvScoringModelFile='" + snvscoringmodelfile + "'") else ""} \
         ~{if defined(indelscoringmodelfile) then ("--indelScoringModelFile='" + indelscoringmodelfile + "'") else ""} \
         ;~{select_first([rundir, "generated"])}/runWorkflow.py \
         ~{if defined(select_first([mode, "local"])) then ("--mode " + select_first([mode, "local"])) else ''} \
         ~{if defined(queue) then ("--queue " + queue) else ''} \
         ~{if defined(memGb) then ("--memGb " + memGb) else ''} \
         ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
         --jobs ~{select_first([runtime_cpu, 4])}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "michaelfranklin/strelka:2.9.10"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4, 4])}G"
       preemptible: 2
     }
     output {
       File configPickle = (select_first([rundir, "generated"]) + "/runWorkflow.py.config.pickle")
       File script = (select_first([rundir, "generated"]) + "/runWorkflow.py")
       File stats = (select_first([rundir, "generated"]) + "/results/stats/runStats.tsv")
       File indels = (select_first([rundir, "generated"]) + "/results/variants/somatic.indels.vcf.gz")
       File indels_tbi = (select_first([rundir, "generated"]) + "/results/variants/somatic.indels.vcf.gz") + ".tbi"
       File snvs = (select_first([rundir, "generated"]) + "/results/variants/somatic.snvs.vcf.gz")
       File snvs_tbi = (select_first([rundir, "generated"]) + "/results/variants/somatic.snvs.vcf.gz") + ".tbi"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Strelka (Somatic)
   doc: |-
     Usage: configureStrelkaSomaticWorkflow.py [options]
     Version: 2.9.10
     This script configures Strelka somatic small variant calling.
     You must specify an alignment file (BAM or CRAM) for each sample of a matched tumor-normal pair.
     Configuration will produce a workflow run script which can execute the workflow on a single node or through
     sge and resume any interrupted execution.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: michaelfranklin/strelka:2.9.10

   inputs:
   - id: normalBam
     label: normalBam
     doc: Normal sample BAM or CRAM file. (no default)
     type: File
     secondaryFiles:
     - pattern: .crai
     inputBinding:
       prefix: --normalBam=
       position: 1
       separate: false
   - id: tumorBam
     label: tumorBam
     doc: (--tumorBam)  Tumor sample BAM or CRAM file. [required] (no default)
     type: File
     secondaryFiles:
     - pattern: .crai
     inputBinding:
       prefix: --tumourBam=
       position: 1
       separate: false
   - id: reference
     label: reference
     doc: ' samtools-indexed reference fasta file [required]'
     type: File
     secondaryFiles:
     - pattern: .fai
     inputBinding:
       prefix: --referenceFasta=
       position: 1
       separate: false
   - id: rundir
     label: rundir
     doc: |-
       Name of directory to be created where all workflow scripts and output will be written. Each analysis requires a separate directory. (default: StrelkaSomaticWorkflow)
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --runDir=
       position: 1
       separate: false
   - id: region
     label: region
     doc: |-
       Limit the analysis to one or more genome region(s) for debugging purposes. If this argument is provided multiple times the union of all specified regions will be analyzed. All regions must be non-overlapping to get a meaningful result. Examples: '--region chr20' (whole chromosome), '--region chr2:100-2000 --region chr3:2500-3000' (two regions)'. If this option is specified (one or more times) together with the 'callRegions' BED file,then all region arguments will be intersected with the callRegions BED track.
     type:
     - type: array
       inputBinding:
         prefix: --region
       items: string
     - 'null'
     inputBinding:
       position: 1
   - id: config
     label: config
     doc: |-
       provide a configuration file to override defaults in global config file (/opt/strelka/bin/configureStrelkaSomaticWorkflow.py.ini)
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --config=
       position: 1
       separate: false
   - id: outputcallableregions
     label: outputcallableregions
     doc: Output a bed file describing somatic callable regions of the genome
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --outputCallableRegions
       position: 1
       separate: true
   - id: indelCandidates
     label: indelCandidates
     doc: |-
       Specify a VCF of candidate indel alleles. These alleles are always evaluated but only reported in the output when they are inferred to exist in the sample. The VCF must be tabix indexed. All indel alleles must be left-shifted/normalized, any unnormalized alleles will be ignored. This option may be specified more than once, multiple input VCFs will be merged. (default: None)
     type:
     - type: array
       inputBinding:
         prefix: --indelCandidates=
         separate: false
       items: File
     - 'null'
     inputBinding:
       position: 1
   - id: forcedgt
     label: forcedgt
     doc: |-
       Specify a VCF of candidate alleles. These alleles are always evaluated and reported even if they are unlikely to exist in the sample. The VCF must be tabix indexed. All indel alleles must be left- shifted/normalized, any unnormalized allele will trigger a runtime error. This option may be specified more than once, multiple input VCFs will be merged. Note that for any SNVs provided in the VCF, the SNV site will be reported (and for gVCF, excluded from block compression), but the specific SNV alleles are ignored. (default: None)
     type:
     - type: array
       inputBinding:
         prefix: --forcedGT=
         separate: false
       items: File
     - 'null'
     inputBinding:
       position: 1
   - id: targeted
     label: targeted
     doc: |-
       Set options for other targeted input: note in particular that this flag turns off high-depth filters
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --targeted
       position: 1
       separate: true
   - id: exome
     label: exome
     doc: |-
       Set options for exome: note in particular that this flag turns off high-depth filters
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --exome
       position: 1
       separate: true
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
       prefix: --callRegions=
       position: 1
       separate: false
   - id: noisevcf
     label: noisevcf
     doc: Noise vcf file (submit argument multiple times for more than one file)
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       prefix: --noiseVcf=
       position: 1
       separate: false
   - id: scansizemb
     label: scansizemb
     doc: |-
       Maximum sequence region size (in megabases) scanned by each task during genome variant calling. (default: 12)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --scanSizeMb=
       position: 1
       separate: false
   - id: callmemmb
     label: callmemmb
     doc: |-
       Set variant calling task memory limit (in megabytes). It is not recommended to change the default in most cases, but this might be required for a sample of unusual depth.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --callMemMb=
       position: 1
       separate: false
   - id: retaintempfiles
     label: retaintempfiles
     doc: Keep all temporary files (for workflow debugging)
     type: boolean
     default: false
     inputBinding:
       prefix: --retainTempFiles
       position: 1
       separate: true
   - id: disableevs
     label: disableevs
     doc: Disable empirical variant scoring (EVS).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --disableEVS
       position: 1
       separate: true
   - id: reportevsfeatures
     label: reportevsfeatures
     doc: ' Report all empirical variant scoring features in VCF output.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --reportEVSFeatures
       position: 1
       separate: true
   - id: snvscoringmodelfile
     label: snvscoringmodelfile
     doc: |2-
        Provide a custom empirical scoring model file for SNVs (default: /opt/strelka/share/config/somaticSNVScoringM odels.json)
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --snvScoringModelFile=
       position: 1
       separate: false
   - id: indelscoringmodelfile
     label: indelscoringmodelfile
     doc: |2-
        Provide a custom empirical scoring model file for indels (default: /opt/strelka/share/config/somaticInde lScoringModels.json)
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --indelScoringModelFile=
       position: 1
       separate: false
   - id: mode
     label: mode
     doc: (-m MODE)  select run mode (local|sge)
     type: string
     default: local
     inputBinding:
       prefix: --mode
       position: 3
       shellQuote: false
   - id: queue
     label: queue
     doc: (-q QUEUE) specify scheduler queue name
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --queue
       position: 3
       shellQuote: false
   - id: memGb
     label: memGb
     doc: |2-
        (-g MEMGB) gigabytes of memory available to run workflow -- only meaningful in local mode, must be an integer (default: Estimate the total memory for this node for local mode, 'unlimited' for sge mode)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --memGb
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

   outputs:
   - id: configPickle
     label: configPickle
     type: File
     outputBinding:
       glob: $((inputs.rundir + "/runWorkflow.py.config.pickle"))
       outputEval: $((inputs.rundir.basename + "/runWorkflow.py.config.pickle"))
       loadContents: false
   - id: script
     label: script
     type: File
     outputBinding:
       glob: $((inputs.rundir + "/runWorkflow.py"))
       outputEval: $((inputs.rundir.basename + "/runWorkflow.py"))
       loadContents: false
   - id: stats
     label: stats
     doc: |-
       A tab-delimited report of various internal statistics from the variant calling process: Runtime information accumulated for each genome segment, excluding auxiliary steps such as BAM indexing and vcf merging. Indel candidacy statistics
     type: File
     outputBinding:
       glob: $((inputs.rundir + "/results/stats/runStats.tsv"))
       outputEval: $((inputs.rundir.basename + "/results/stats/runStats.tsv"))
       loadContents: false
   - id: indels
     label: indels
     doc: ''
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $((inputs.rundir + "/results/variants/somatic.indels.vcf.gz"))
       outputEval: $((inputs.rundir.basename + "/results/variants/somatic.indels.vcf.gz"))
       loadContents: false
   - id: snvs
     label: snvs
     doc: ''
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $((inputs.rundir + "/results/variants/somatic.snvs.vcf.gz"))
       outputEval: $((inputs.rundir.basename + "/results/variants/somatic.snvs.vcf.gz"))
       loadContents: false
   stdout: _stdout
   stderr: _stderr
   arguments:
   - position: 0
     valueFrom: configureStrelkaSomaticWorkflow.py
   - position: 2
     valueFrom: $(";{rundir}/runWorkflow.py".replace(/\{rundir\}/g, inputs.rundir))
     shellQuote: false
   - prefix: --jobs
     position: 3
     valueFrom: $([inputs.runtime_cpu, 4].filter(function (inner) { return inner != null
       })[0])
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: strelka_somatic_cram


