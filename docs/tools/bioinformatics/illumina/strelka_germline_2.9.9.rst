:orphan:

Strelka (Germline)
=====================================

``strelka_germline`` · *1 contributor · 2 versions*

Strelka2 is a fast and accurate small variant caller optimized for analysis of germline variation 
in small cohorts and somatic variation in tumor/normal sample pairs. The germline caller employs 
an efficient tiered haplotype model to improve accuracy and provide read-backed phasing, adaptively 
selecting between assembly and a faster alignment-based haplotyping approach at each variant locus. 
The germline caller also analyzes input sequencing data using a mixture-model indel error estimation 
method to improve robustness to indel noise. The somatic calling model improves on the original 
Strelka method for liquid and late-stage tumor analysis by accounting for possible tumor cell 
contamination in the normal sample. A final empirical variant re-scoring step using random forest 
models trained on various call quality features has been added to both callers to further improve precision.

Compared with submissions to the recent PrecisonFDA Consistency and Truth challenges, the average 
indel F-score for Strelka2 running in its default configuration is 3.1% and 0.08% higher, respectively, 
than the best challenge submissions. Runtime on a 28-core server is ~40 minutes for 40x WGS germline 
analysis and ~3 hours for a 110x/40x WGS tumor-normal somatic analysis

Strelka accepts input read mappings from BAM or CRAM files, and optionally candidate and/or forced-call 
alleles from VCF. It reports all small variant predictions in VCF 4.1 format. Germline variant 
reporting uses the gVCF conventions to represent both variant and reference call confidence. 
For best somatic indel performance, Strelka is designed to be run with the Manta structural variant 
and indel caller, which provides additional indel candidates up to a given maxiumum indel size 
(49 by default). By design, Manta and Strelka run together with default settings provide complete 
coverage over all indel sizes (in additional to SVs and SNVs). 

See the user guide for a full description of capabilities and limitations

.. warning::

   Strelka (Germline) did not include a container. You can provide one through the command line by including
   the following instruction:

   .. code-block:: bash

      janis run --container-override 'strelka_germline=<organisation/container:version>' strelka_germline
    
Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.illumina.strelkagermline.strelkagermline import StrelkaGermline_2_9_9

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "strelka_germline_step",
           StrelkaGermline_2_9_9(
               bam=None,
               reference=None,
           )
       )
       wf.output("configPickle", source=strelka_germline_step.configPickle)
       wf.output("script", source=strelka_germline_step.script)
       wf.output("stats", source=strelka_germline_step.stats)
       wf.output("variants", source=strelka_germline_step.variants)
       wf.output("genome", source=strelka_germline_step.genome)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for strelka_germline:

.. code-block:: bash

   # user inputs
   janis inputs strelka_germline > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       reference: reference.fasta




5. Run strelka_germline with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       --container-override 'strelka_germline=<organisation/container:version>' \
       strelka_germline





Information
------------

:ID: ``strelka_germline``
:URL: `https://github.com/Illumina/strelka <https://github.com/Illumina/strelka>`_
:Versions: 2.9.10, 2.9.9
:Container: 
:Authors: Michael Franklin
:Citations: None
:Created: 2018-12-24
:Updated: 2019-01-24


Outputs
-----------

============  ============  ===========================================================================================================================================================================================================================================
name          type          documentation
============  ============  ===========================================================================================================================================================================================================================================
configPickle  File
script        File
stats         tsv           A tab-delimited report of various internal statistics from the variant calling process: Runtime information accumulated for each genome segment, excluding auxiliary steps such as BAM indexing and vcf merging. Indel candidacy statistics
variants      Gzipped<VCF>  Primary variant inferences are provided as a series of VCF 4.1 files
genome        Gzipped<VCF>
============  ============  ===========================================================================================================================================================================================================================================


Additional configuration (inputs)
---------------------------------

========================  ======================  ==================  ==========  ====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                      type                    prefix                position  documentation
========================  ======================  ==================  ==========  ====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
bam                       IndexedBam              --bam                        1  Sample BAM or CRAM file. May be specified more than once, multiple inputs will be treated as each BAM file representing a different sample. [required] (no default)
reference                 FastaWithIndexes        --referenceFasta             1  samtools-indexed reference fasta file [required]
relativeStrelkaDirectory  Optional<String>        --runDir                     1  Name of directory to be created where all workflow scripts and output will be written. Each analysis requires a separate directory.
ploidy                    Optional<Gzipped<VCF>>  --ploidy                     1  Provide ploidy file in VCF. The VCF should include one sample column per input sample labeled with the same sample names found in the input BAM/CRAM RG header sections. Ploidy should be provided in records using the FORMAT/CN field, which are interpreted to span the range [POS+1, INFO/END]. Any CN value besides 1 or 0 will be treated as 2. File must be tabix indexed. (no default)
noCompress                Optional<Gzipped<VCF>>  --noCompress                 1  Provide BED file of regions where gVCF block compression is not allowed. File must be bgzip- compressed/tabix-indexed. (no default)
callContinuousVf          Optional<String>        --callContinuousVf              Call variants on CHROM without a ploidy prior assumption, issuing calls with continuous variant frequencies (no default)
rna                       Optional<Boolean>       --rna                        1  Set options for RNA-Seq input.
indelCandidates           Optional<Gzipped<VCF>>  --indelCandidates            1  Specify a VCF of candidate indel alleles. These alleles are always evaluated but only reported in the output when they are inferred to exist in the sample. The VCF must be tabix indexed. All indel alleles must be left-shifted/normalized, any unnormalized alleles will be ignored. This option may be specified more than once, multiple input VCFs will be merged. (default: None)
forcedGT                  Optional<Gzipped<VCF>>  --forcedGT                   1  Specify a VCF of candidate alleles. These alleles are always evaluated and reported even if they are unlikely to exist in the sample. The VCF must be tabix indexed. All indel alleles must be left- shifted/normalized, any unnormalized allele will trigger a runtime error. This option may be specified more than once, multiple input VCFs will be merged. Note that for any SNVs provided in the VCF, the SNV site will be reported (and for gVCF, excluded from block compression), but the specific SNV alleles are ignored. (default: None)
exome                     Optional<Boolean>       --exome                      1  Set options for exome note in particular that this flag turns off high-depth filters
targeted                  Optional<Boolean>       --exome                      1  Set options for other targeted input: note in particular that this flag turns off high-depth filters
callRegions               Optional<Gzipped<bed>>  --callRegions=               1  Optionally provide a bgzip-compressed/tabix-indexed BED file containing the set of regions to call. No VCF output will be provided outside of these regions. The full genome will still be used to estimate statistics from the input (such as expected depth per chromosome). Only one BED file may be specified. (default: call the entire genome)
mode                      Optional<String>        --mode                       3  (-m MODE)  select run mode (local|sge)
queue                     Optional<String>        --queue                      3  (-q QUEUE) specify scheduler queue name
memGb                     Optional<String>        --memGb                      3  (-g MEMGB) gigabytes of memory available to run workflow -- only meaningful in local mode, must be an integer (default: Estimate the total memory for this node for local mode, 'unlimited' for sge mode)
quiet                     Optional<Boolean>       --quiet                      3  Don't write any log output to stderr (but still write to workspace/pyflow.data/logs/pyflow_log.txt)
mailTo                    Optional<String>        --mailTo                     3  (-e) send email notification of job completion status to this address (may be provided multiple times for more than one email address)
========================  ======================  ==================  ==========  ====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task strelka_germline {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File bam
       File bam_bai
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       String? relativeStrelkaDirectory
       File? ploidy
       File? ploidy_tbi
       File? noCompress
       File? noCompress_tbi
       String? callContinuousVf
       Boolean? rna
       File? indelCandidates
       File? indelCandidates_tbi
       File? forcedGT
       File? forcedGT_tbi
       Boolean? exome
       Boolean? targeted
       File? callRegions
       File? callRegions_tbi
       String? mode
       String? queue
       String? memGb
       Boolean? quiet
       String? mailTo
     }
     command <<<
       set -e
        \
         ~{if defined(callContinuousVf) then ("--callContinuousVf '" + callContinuousVf + "'") else ""} \
         configureStrelkaGermlineWorkflow.py \
         --bam ~{bam} \
         --referenceFasta ~{reference} \
         ~{if defined(select_first([relativeStrelkaDirectory, "strelka_dir"])) then ("--runDir " + select_first([relativeStrelkaDirectory, "strelka_dir"])) else ''} \
         ~{if defined(ploidy) then ("--ploidy " + ploidy) else ''} \
         ~{if defined(noCompress) then ("--noCompress " + noCompress) else ''} \
         ~{if (defined(rna) && select_first([rna])) then "--rna" else ""} \
         ~{if defined(indelCandidates) then ("--indelCandidates " + indelCandidates) else ''} \
         ~{if defined(forcedGT) then ("--forcedGT " + forcedGT) else ''} \
         ~{if (defined(exome) && select_first([exome])) then "--exome" else ""} \
         ~{if (defined(targeted) && select_first([targeted])) then "--exome" else ""} \
         ~{if defined(callRegions) then ("--callRegions='" + callRegions + "'") else ""} \
         ;~{select_first([relativeStrelkaDirectory, "strelka_dir"])}/runWorkflow.py \
         ~{if defined(select_first([mode, "local"])) then ("--mode " + select_first([mode, "local"])) else ''} \
         ~{if defined(queue) then ("--queue " + queue) else ''} \
         ~{if defined(memGb) then ("--memGb " + memGb) else ''} \
         ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
         ~{if defined(mailTo) then ("--mailTo " + mailTo) else ''} \
         --jobs ~{select_first([runtime_cpu, 4])}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: ""
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4, 4])}G"
       preemptible: 2
     }
     output {
       File configPickle = (select_first([relativeStrelkaDirectory, "strelka_dir"]) + "/runWorkflow.py.config.pickle")
       File script = (select_first([relativeStrelkaDirectory, "strelka_dir"]) + "/runWorkflow.py")
       File stats = (select_first([relativeStrelkaDirectory, "strelka_dir"]) + "/results/stats/runStats.tsv")
       File variants = (select_first([relativeStrelkaDirectory, "strelka_dir"]) + "/results/variants/variants.vcf.gz")
       File variants_tbi = (select_first([relativeStrelkaDirectory, "strelka_dir"]) + "/results/variants/variants.vcf.gz") + ".tbi"
       File genome = (select_first([relativeStrelkaDirectory, "strelka_dir"]) + "/results/variants/genome.vcf.gz")
       File genome_tbi = (select_first([relativeStrelkaDirectory, "strelka_dir"]) + "/results/variants/genome.vcf.gz") + ".tbi"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Strelka (Germline)
   doc: |-
     Strelka2 is a fast and accurate small variant caller optimized for analysis of germline variation 
     in small cohorts and somatic variation in tumor/normal sample pairs. The germline caller employs 
     an efficient tiered haplotype model to improve accuracy and provide read-backed phasing, adaptively 
     selecting between assembly and a faster alignment-based haplotyping approach at each variant locus. 
     The germline caller also analyzes input sequencing data using a mixture-model indel error estimation 
     method to improve robustness to indel noise. The somatic calling model improves on the original 
     Strelka method for liquid and late-stage tumor analysis by accounting for possible tumor cell 
     contamination in the normal sample. A final empirical variant re-scoring step using random forest 
     models trained on various call quality features has been added to both callers to further improve precision.

     Compared with submissions to the recent PrecisonFDA Consistency and Truth challenges, the average 
     indel F-score for Strelka2 running in its default configuration is 3.1% and 0.08% higher, respectively, 
     than the best challenge submissions. Runtime on a 28-core server is ~40 minutes for 40x WGS germline 
     analysis and ~3 hours for a 110x/40x WGS tumor-normal somatic analysis

     Strelka accepts input read mappings from BAM or CRAM files, and optionally candidate and/or forced-call 
     alleles from VCF. It reports all small variant predictions in VCF 4.1 format. Germline variant 
     reporting uses the gVCF conventions to represent both variant and reference call confidence. 
     For best somatic indel performance, Strelka is designed to be run with the Manta structural variant 
     and indel caller, which provides additional indel candidates up to a given maxiumum indel size 
     (49 by default). By design, Manta and Strelka run together with default settings provide complete 
     coverage over all indel sizes (in additional to SVs and SNVs). 

     See the user guide for a full description of capabilities and limitations

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ''

   inputs:
   - id: bam
     label: bam
     doc: |-
       Sample BAM or CRAM file. May be specified more than once, multiple inputs will be treated as each BAM file representing a different sample. [required] (no default)
     type: File
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       prefix: --bam
       position: 1
       shellQuote: false
   - id: reference
     label: reference
     doc: samtools-indexed reference fasta file [required]
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
       prefix: --referenceFasta
       position: 1
       shellQuote: false
   - id: relativeStrelkaDirectory
     label: relativeStrelkaDirectory
     doc: |-
       Name of directory to be created where all workflow scripts and output will be written. Each analysis requires a separate directory.
     type: string
     default: strelka_dir
     inputBinding:
       prefix: --runDir
       position: 1
       shellQuote: false
   - id: ploidy
     label: ploidy
     doc: |-
       Provide ploidy file in VCF. The VCF should include one sample column per input sample labeled with the same sample names found in the input BAM/CRAM RG header sections. Ploidy should be provided in records using the FORMAT/CN field, which are interpreted to span the range [POS+1, INFO/END]. Any CN value besides 1 or 0 will be treated as 2. File must be tabix indexed. (no default)
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       prefix: --ploidy
       position: 1
       shellQuote: false
   - id: noCompress
     label: noCompress
     doc: |-
       Provide BED file of regions where gVCF block compression is not allowed. File must be bgzip- compressed/tabix-indexed. (no default)
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       prefix: --noCompress
       position: 1
       shellQuote: false
   - id: callContinuousVf
     label: callContinuousVf
     doc: |-
       Call variants on CHROM without a ploidy prior assumption, issuing calls with continuous variant frequencies (no default)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --callContinuousVf
   - id: rna
     label: rna
     doc: Set options for RNA-Seq input.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --rna
       position: 1
       shellQuote: false
   - id: indelCandidates
     label: indelCandidates
     doc: |-
       Specify a VCF of candidate indel alleles. These alleles are always evaluated but only reported in the output when they are inferred to exist in the sample. The VCF must be tabix indexed. All indel alleles must be left-shifted/normalized, any unnormalized alleles will be ignored. This option may be specified more than once, multiple input VCFs will be merged. (default: None)
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       prefix: --indelCandidates
       position: 1
       shellQuote: false
   - id: forcedGT
     label: forcedGT
     doc: |-
       Specify a VCF of candidate alleles. These alleles are always evaluated and reported even if they are unlikely to exist in the sample. The VCF must be tabix indexed. All indel alleles must be left- shifted/normalized, any unnormalized allele will trigger a runtime error. This option may be specified more than once, multiple input VCFs will be merged. Note that for any SNVs provided in the VCF, the SNV site will be reported (and for gVCF, excluded from block compression), but the specific SNV alleles are ignored. (default: None)
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       prefix: --forcedGT
       position: 1
       shellQuote: false
   - id: exome
     label: exome
     doc: |-
       Set options for exome note in particular that this flag turns off high-depth filters
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --exome
       position: 1
       shellQuote: false
   - id: targeted
     label: targeted
     doc: |-
       Set options for other targeted input: note in particular that this flag turns off high-depth filters
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --exome
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
       prefix: --callRegions=
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
   - id: mailTo
     label: mailTo
     doc: |-
       (-e) send email notification of job completion status to this address (may be provided multiple times for more than one email address)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --mailTo
       position: 3
       shellQuote: false

   outputs:
   - id: configPickle
     label: configPickle
     type: File
     outputBinding:
       glob: $((inputs.relativeStrelkaDirectory + "/runWorkflow.py.config.pickle"))
       outputEval: $((inputs.relativeStrelkaDirectory.basename + "/runWorkflow.py.config.pickle"))
       loadContents: false
   - id: script
     label: script
     type: File
     outputBinding:
       glob: $((inputs.relativeStrelkaDirectory + "/runWorkflow.py"))
       outputEval: $((inputs.relativeStrelkaDirectory.basename + "/runWorkflow.py"))
       loadContents: false
   - id: stats
     label: stats
     doc: |-
       A tab-delimited report of various internal statistics from the variant calling process: Runtime information accumulated for each genome segment, excluding auxiliary steps such as BAM indexing and vcf merging. Indel candidacy statistics
     type: File
     outputBinding:
       glob: $((inputs.relativeStrelkaDirectory + "/results/stats/runStats.tsv"))
       outputEval: $((inputs.relativeStrelkaDirectory.basename + "/results/stats/runStats.tsv"))
       loadContents: false
   - id: variants
     label: variants
     doc: Primary variant inferences are provided as a series of VCF 4.1 files
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $((inputs.relativeStrelkaDirectory + "/results/variants/variants.vcf.gz"))
       outputEval: |-
         $((inputs.relativeStrelkaDirectory.basename + "/results/variants/variants.vcf.gz"))
       loadContents: false
   - id: genome
     label: genome
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $((inputs.relativeStrelkaDirectory + "/results/variants/genome.vcf.gz"))
       outputEval: |-
         $((inputs.relativeStrelkaDirectory.basename + "/results/variants/genome.vcf.gz"))
       loadContents: false
   stdout: _stdout
   stderr: _stderr
   arguments:
   - position: 0
     valueFrom: configureStrelkaGermlineWorkflow.py
     shellQuote: false
   - position: 2
     valueFrom: |-
       $(";{relativeStrelkaDirectory}/runWorkflow.py".replace(/\{relativeStrelkaDirectory\}/g, inputs.relativeStrelkaDirectory))
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
   id: strelka_germline


