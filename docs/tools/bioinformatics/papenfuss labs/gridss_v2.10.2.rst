:orphan:

Gridss
===============

``gridss`` · *1 contributor · 7 versions*

GRIDSS: the Genomic Rearrangement IDentification Software Suite

GRIDSS is a module software suite containing tools useful for the detection of genomic rearrangements.
GRIDSS includes a genome-wide break-end assembler, as well as a structural variation caller for Illumina
sequencing data. GRIDSS calls variants based on alignment-guided positional de Bruijn graph genome-wide
break-end assembly, split read, and read pair evidence.

GRIDSS makes extensive use of the standard tags defined by SAM specifications. Due to the modular design,
any step (such as split read identification) can be replaced by another implementation that also outputs
using the standard tags. It is hoped that GRIDSS can serve as an exemplar modular structural variant
pipeline designed for interoperability with other tools.

If you have any trouble running GRIDSS, please raise an issue using the Issues tab above. Based on feedback
from users, a user guide will be produced outlining common workflows, pitfalls, and use cases.



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.papenfuss.gridss.gridss import Gridss_2_10_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gridss_step",
           Gridss_2_10_2(
               bams=None,
               reference=None,
           )
       )
       wf.output("vcf", source=gridss_step.vcf)
       wf.output("assembly", source=gridss_step.assembly)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for gridss:

.. code-block:: bash

   # user inputs
   janis inputs gridss > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bams:
       - bams_0.bam
       - bams_1.bam
       reference: reference.fasta




5. Run gridss with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       gridss

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          gridss

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``gridss``
:URL: `https://github.com/PapenfussLab/gridss/wiki/GRIDSS-Documentation <https://github.com/PapenfussLab/gridss/wiki/GRIDSS-Documentation>`_
:Versions: v2.9.4, v2.8.3, v2.6.2, v2.5.1-dev, v2.4.0, v2.2.3, v2.10.2
:Container: gridss/gridss:2.10.2
:Authors: Jiaan Yu
:Citations: Daniel L. Cameron, Jan Schröder, Jocelyn Sietsma Penington, Hongdo Do, Ramyar Molania, Alexander Dobrovic, Terence P. Speed and Anthony T. Papenfuss. GRIDSS: sensitive and specific genomic rearrangement detection using positional de Bruijn graph assembly. Genome Research, 2017 doi: 10.1101/gr.222109.117
:DOI: 10.1101/gr.222109.117
:Created: 2021-03-30
:Updated: 2021-03-30


Outputs
-----------

========  ======  ===============
name      type    documentation
========  ======  ===============
vcf       VCF
assembly  BAM
========  ======  ===============


Additional configuration (inputs)
---------------------------------

==============================  =======================  ================================  ==========  ==================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                            type                     prefix                              position  documentation
==============================  =======================  ================================  ==========  ==================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
bams                            Array<IndexedBam>                                                  10
reference                       FastaWithIndexes         --reference                                   reference genome to use.
outputFilename                  Optional<Filename>       --output                                      output VCF.
assemblyFilename                Optional<Filename>       --assembly                                    location of the GRIDSS assembly BAM. This file will be created by GRIDSS.
threads                         Optional<Integer>        --threads                                     number of threads to use. (Default: 8)
jarPath                         Optional<String>         --jar                                         location of GRIDSS jar
workingDir                      Optional<String>         --workingdir                                  directory to place GRIDSS intermediate and temporary files. .gridss.working subdirectories will be created. (Default: .)
blacklist                       Optional<bed>            --blacklist                                   BED file containing regions to ignore
steps                           Optional<Array<String>>  --steps                                       processing steps to run. Defaults to all steps. Multiple steps are specified using comma separators. Possible steps are: setupreference, preprocess, assemble, call, all. WARNING: multiple instances of GRIDSS generating reference files at the same time will result in file corruption. Make sure these files are generated before runninng parallel GRIDSS jobs.
configuration                   Optional<File>           --configuration                               configuration file use to override default GRIDSS settings.
labels                          Optional<Array<String>>  --labels                                      comma separated labels to use in the output VCF for the input files. Supporting read counts for input files with the same label are aggregated (useful for multiple sequencing runs of the same sample). Labels default to input filenames, unless a single read group with a non-empty sample name exists in which case the read group sample name is used (which can be disabled by "useReadGroupSampleNameCategoryLabel=false" in the configuration file). If labels are specified, they must be specified for all input files.
externalaligner                 Optional<String>         --externalaligner                             use the system version of bwa instead of the in-process version packaged with GRIDSS
jvmheap                         Optional<String>         --jvmheap                                     size of JVM heap for assembly and variant calling. (Default: 30g)
maxcoverage                     Optional<Integer>        --maxcoverage                                 maximum coverage. Regions with coverage in excess of this are ignored. (Default: 50000)
picardoptions                   Optional<String>         --picardoptions                               additional standard Picard command line options. Useful options include VALIDATION_STRINGENCY=LENIENT and COMPRESSION_LEVEL=0. See https://broadinstitute.github.io/picard/command-line-overview.html
useproperpair                   Optional<String>         --useproperpair                               use SAM 'proper pair' flag to determine whether a read pair is discordant. Default: use library fragment size distribution to determine read pair concordance
concordantreadpairdistribution  Optional<Float>          --concordantreadpairdistribution              portion of 6 sigma read pairs distribution considered concordantly mapped. (Default: 0.995)
keepTempFiles                   Optional<Boolean>        --keepTempFiles                               keep intermediate files. Not recommended except for debugging due to the high disk usage.
nojni                           Optional<Boolean>        --nojni                                       do not use JNI native code acceleration libraries (snappy, GKL, ssw, bwa).
jobindex                        Optional<Integer>        --jobindex                                    zero-based assembly job index (only required when performing parallel assembly across multiple computers)
jobnodes                        Optional<Integer>        --jobnodes                                    total number of assembly jobs (only required when performing parallel assembly across multiple computers). Note than an assembly job with any --job argument is required to be run after all indexed jobs have been completed to gather the output files together.
==============================  =======================  ================================  ==========  ==================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task gridss {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       Array[File] bams
       Array[File] bams_bai
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       String? outputFilename
       String? assemblyFilename
       Int? threads
       String? jarPath
       String? workingDir
       File? blacklist
       Array[String]? steps
       File? configuration
       Array[String]? labels
       String? externalaligner
       String? jvmheap
       Int? maxcoverage
       String? picardoptions
       String? useproperpair
       Float? concordantreadpairdistribution
       Boolean? keepTempFiles
       Boolean? nojni
       Int? jobindex
       Int? jobnodes
     }

     command <<<
       set -e
       /opt/gridss/gridss.sh \
         --reference '~{reference}' \
         --output '~{select_first([outputFilename, "generated.vcf"])}' \
         --assembly '~{select_first([assemblyFilename, "generated.assembly.bam"])}' \
         ~{if defined(select_first([threads, select_first([runtime_cpu, 1])])) then ("--threads " + select_first([threads, select_first([runtime_cpu, 1])])) else ''} \
         ~{if defined(jarPath) then ("--jar '" + jarPath + "'") else ""} \
         ~{if defined(select_first([workingDir, "./TMP"])) then ("--workingdir '" + select_first([workingDir, "./TMP"]) + "'") else ""} \
         ~{if defined(blacklist) then ("--blacklist '" + blacklist + "'") else ""} \
         ~{if (defined(steps) && length(select_first([steps])) > 0) then "--steps '" + sep("','", select_first([steps])) + "'" else ""} \
         ~{if defined(configuration) then ("--configuration '" + configuration + "'") else ""} \
         ~{if (defined(labels) && length(select_first([labels])) > 0) then "--labels '" + sep("','", select_first([labels])) + "'" else ""} \
         ~{if defined(externalaligner) then ("--externalaligner '" + externalaligner + "'") else ""} \
         ~{if defined(jvmheap) then ("--jvmheap '" + jvmheap + "'") else ""} \
         ~{if defined(maxcoverage) then ("--maxcoverage " + maxcoverage) else ''} \
         ~{if defined(picardoptions) then ("--picardoptions '" + picardoptions + "'") else ""} \
         ~{if defined(useproperpair) then ("--useproperpair '" + useproperpair + "'") else ""} \
         ~{if defined(concordantreadpairdistribution) then ("--concordantreadpairdistribution " + concordantreadpairdistribution) else ''} \
         ~{if (defined(keepTempFiles) && select_first([keepTempFiles])) then "--keepTempFiles" else ""} \
         ~{if (defined(nojni) && select_first([nojni])) then "--nojni" else ""} \
         ~{if defined(jobindex) then ("--jobindex " + jobindex) else ''} \
         ~{if defined(jobnodes) then ("--jobnodes " + jobnodes) else ''} \
         ~{if length(bams) > 0 then "'" + sep("' '", bams) + "'" else ""}
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 8, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "gridss/gridss:2.10.2"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 31, 4])}G"
       preemptible: 2
     }

     output {
       File vcf = select_first([outputFilename, "generated.vcf"])
       File assembly = select_first([assemblyFilename, "generated.assembly.bam"])
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Gridss

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: gridss/gridss:2.10.2

   inputs:
   - id: bams
     label: bams
     type:
       type: array
       items: File
     inputBinding:
       position: 10
   - id: reference
     label: reference
     doc: reference genome to use.
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
       prefix: --reference
   - id: outputFilename
     label: outputFilename
     doc: output VCF.
     type:
     - string
     - 'null'
     default: generated.vcf
     inputBinding:
       prefix: --output
   - id: assemblyFilename
     label: assemblyFilename
     doc: location of the GRIDSS assembly BAM. This file will be created by GRIDSS.
     type:
     - string
     - 'null'
     default: generated.assembly.bam
     inputBinding:
       prefix: --assembly
   - id: threads
     label: threads
     doc: 'number of threads to use. (Default: 8)'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --threads
       valueFrom: |-
         $([inputs.runtime_cpu, 8, 1].filter(function (inner) { return inner != null })[0])
   - id: jarPath
     label: jarPath
     doc: location of GRIDSS jar
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --jar
   - id: workingDir
     label: workingDir
     doc: |-
       directory to place GRIDSS intermediate and temporary files. .gridss.working subdirectories will be created. (Default: .)
     type: string
     default: ./TMP
     inputBinding:
       prefix: --workingdir
   - id: blacklist
     label: blacklist
     doc: BED file containing regions to ignore
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --blacklist
   - id: steps
     label: steps
     doc: |-
       processing steps to run. Defaults to all steps. Multiple steps are specified using comma separators. Possible steps are: setupreference, preprocess, assemble, call, all. WARNING: multiple instances of GRIDSS generating reference files at the same time will result in file corruption. Make sure these files are generated before runninng parallel GRIDSS jobs.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --steps
       itemSeparator: ','
   - id: configuration
     label: configuration
     doc: configuration file use to override default GRIDSS settings.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --configuration
   - id: labels
     label: labels
     doc: |-
       comma separated labels to use in the output VCF for the input files. Supporting read counts for input files with the same label are aggregated (useful for multiple sequencing runs of the same sample). Labels default to input filenames, unless a single read group with a non-empty sample name exists in which case the read group sample name is used (which can be disabled by "useReadGroupSampleNameCategoryLabel=false" in the configuration file). If labels are specified, they must be specified for all input files.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --labels
       itemSeparator: ','
   - id: externalaligner
     label: externalaligner
     doc: |-
       use the system version of bwa instead of the in-process version packaged with GRIDSS
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --externalaligner
   - id: jvmheap
     label: jvmheap
     doc: 'size of JVM heap for assembly and variant calling. (Default: 30g)'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --jvmheap
   - id: maxcoverage
     label: maxcoverage
     doc: |-
       maximum coverage. Regions with coverage in excess of this are ignored. (Default: 50000)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --maxcoverage
   - id: picardoptions
     label: picardoptions
     doc: |-
       additional standard Picard command line options. Useful options include VALIDATION_STRINGENCY=LENIENT and COMPRESSION_LEVEL=0. See https://broadinstitute.github.io/picard/command-line-overview.html
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --picardoptions
   - id: useproperpair
     label: useproperpair
     doc: |-
       use SAM 'proper pair' flag to determine whether a read pair is discordant. Default: use library fragment size distribution to determine read pair concordance
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --useproperpair
   - id: concordantreadpairdistribution
     label: concordantreadpairdistribution
     doc: |-
       portion of 6 sigma read pairs distribution considered concordantly mapped. (Default: 0.995)
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --concordantreadpairdistribution
   - id: keepTempFiles
     label: keepTempFiles
     doc: |-
       keep intermediate files. Not recommended except for debugging due to the high disk usage.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --keepTempFiles
   - id: nojni
     label: nojni
     doc: do not use JNI native code acceleration libraries (snappy, GKL, ssw, bwa).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --nojni
   - id: jobindex
     label: jobindex
     doc: |-
       zero-based assembly job index (only required when performing parallel assembly across multiple computers)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --jobindex
   - id: jobnodes
     label: jobnodes
     doc: |-
       total number of assembly jobs (only required when performing parallel assembly across multiple computers). Note than an assembly job with any --job argument is required to be run after all indexed jobs have been completed to gather the output files together.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --jobnodes

   outputs:
   - id: vcf
     label: vcf
     type: File
     outputBinding:
       glob: generated.vcf
       loadContents: false
   - id: assembly
     label: assembly
     type: File
     outputBinding:
       glob: generated.assembly.bam
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: /opt/gridss/gridss.sh
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: gridss


