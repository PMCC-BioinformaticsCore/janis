:orphan:

Gridss
===============

``gridss`` · *1 contributor · 5 versions*

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

       from janis_bioinformatics.tools.papenfuss.gridss.gridss import Gridss_2_2_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gridss_step",
           Gridss_2_2_3(
               reference=None,
               bams=None,
           )
       )
       wf.output("vcf", source=gridss_step.vcf)
       wf.output("assembly", source=gridss_step.assembly)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



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





Information
------------

:ID: ``gridss``
:URL: `https://github.com/PapenfussLab/gridss/wiki/GRIDSS-Documentation <https://github.com/PapenfussLab/gridss/wiki/GRIDSS-Documentation>`_
:Versions: v2.9.4, v2.6.2, v2.5.1-dev, v2.4.0, v2.2.3
:Container: gridss/gridss:v2.2.3
:Authors: Michael Franklin
:Citations: Daniel L. Cameron, Jan Schröder, Jocelyn Sietsma Penington, Hongdo Do, Ramyar Molania, Alexander Dobrovic, Terence P. Speed and Anthony T. Papenfuss. GRIDSS: sensitive and specific genomic rearrangement detection using positional de Bruijn graph assembly. Genome Research, 2017 doi: 10.1101/gr.222109.117
:DOI: 10.1101/gr.222109.117
:Created: 2019-06-19
:Updated: 2019-07-03


Outputs
-----------

========  ==========  ===============
name      type        documentation
========  ==========  ===============
vcf       VCF
assembly  IndexedBam
========  ==========  ===============


Additional configuration (inputs)
---------------------------------

=========================  ==================  =============================  ==========  ===================================================================================================================================================================================================================================================================================================================================
name                       type                prefix                         position    documentation
=========================  ==================  =============================  ==========  ===================================================================================================================================================================================================================================================================================================================================
reference                  FastaWithIndexes    REFERENCE_SEQUENCE=
bams                       Array<IndexedBam>   INPUT=                                     (I=File Coordinate-sorted input BAM file. Default value: null. This option may be specified 0 or more times.
outputFilename             Optional<Filename>  OUTPUT=                                    (O=) VCF structural variation calls. Required.
assemblyFilename           Optional<Filename>  ASSEMBLY=                                  Breakend assemblies which have undergone split read identification Required.
inputLabel                 Optional<String>    INPUT_LABEL=                               Input label. Variant calling evidence breakdowns are reported for each label. Default labels correspond to INPUT filenames. When specifying labels, labels must be provided for all input files. Default value: null. This option may be specified 0 or more times.
inputMaxFragmentSize       Optional<Integer>   INPUT_MAX_FRAGMENT_SIZE=                   Per input maximum concordant fragment size. Default value: null. This option may be specified 0 or more times.
inputMinFragmentSize       Optional<Integer>   INPUT_MIN_FRAGMENT_SIZE=                   Per input minimum concordant fragment size. Default value: null. This option may be specified 0 or more times.
readPairConcordantPercent  Optional<Float>     READ_PAIR_CONCORDANT_PERCENT=              Percent of read pairs considered concorant (0.0-1.0). If this is unset, the SAM proper pair flag is used to determine whether a read is discordantly aligned. Explicit fragment size specification overrides this setting. Default value: 0.995. This option can be set to 'null' to clear the default value.
blacklist                  Optional<bed>       BLACKLIST=                                 (BL=File) BED blacklist of regions to ignore. Assembly of regions such as high-coverage centromeric repeats is slow, and if such regions are to be filtered in downstream analysis anyway, blacklisting those region will improve runtime performance. For human WGS, the ENCODE DAC blacklist is recommended. Default value: null.
configurationFile          Optional<File>      CONFIGURATION_FILE=                        (C=File) gridss configuration file containing overrides Default value: null.
workerThreads              Optional<Integer>   WORKER_THREADS=                            (THREADS=Integer  Number of worker threads to spawn. Defaults to number of cores available. Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread. Default value: 6. This option can be set to 'null' to clear the default value.
workingDir                 Optional<String>    WORKING_DIR=                               Directory to place intermediate results directories. Default location is the same directory as the associated input or output file. Default value: null.
ignoreDuplicates           Optional<Boolean>   IGNORE_DUPLICATES=                         Ignore reads marked as duplicates. Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false}
=========================  ==================  =============================  ==========  ===================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task gridss {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       String? outputFilename
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       Array[File] bams
       Array[File] bams_bai
       String? assemblyFilename
       String? inputLabel
       Int? inputMaxFragmentSize
       Int? inputMinFragmentSize
       Float? readPairConcordantPercent
       File? blacklist
       File? configurationFile
       Int? workerThreads
       String? workingDir
       Boolean? ignoreDuplicates
     }
     command <<<
       set -e
       java -XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap -XX:MaxRAMFraction=1 -XshowSettings:vm -Dsamjdk.create_index=true -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=true -Dgridss.gridss.output_to_temp_file=true -Dsamjdk.buffer_size=4194304 -cp /data/gridss/gridss-2.2.3-gridss-jar-with-dependencies.jar gridss.CallVariants \
         OUTPUT='~{select_first([outputFilename, "generated.vcf"])}' \
         REFERENCE_SEQUENCE='~{reference}' \
         ~{if length(bams) > 0 then "INPUT='" + sep("' INPUT='", bams) + "'" else ""} \
         ASSEMBLY='~{select_first([assemblyFilename, "generated.assembled.bam"])}' \
         ~{if defined(inputLabel) then ("INPUT_LABEL='" + inputLabel + "'") else ""} \
         ~{if defined(inputMaxFragmentSize) then ("INPUT_MAX_FRAGMENT_SIZE=" + inputMaxFragmentSize) else ''} \
         ~{if defined(inputMinFragmentSize) then ("INPUT_MIN_FRAGMENT_SIZE=" + inputMinFragmentSize) else ''} \
         ~{if defined(readPairConcordantPercent) then ("READ_PAIR_CONCORDANT_PERCENT=" + readPairConcordantPercent) else ''} \
         ~{if defined(blacklist) then ("BLACKLIST='" + blacklist + "'") else ""} \
         ~{if defined(configurationFile) then ("CONFIGURATION_FILE='" + configurationFile + "'") else ""} \
         ~{if defined(workerThreads) then ("WORKER_THREADS=" + workerThreads) else ''} \
         ~{if defined(select_first([workingDir, "."])) then ("WORKING_DIR='" + select_first([workingDir, "."]) + "'") else ""} \
         ~{if (defined(ignoreDuplicates) && select_first([ignoreDuplicates])) then "IGNORE_DUPLICATES=" else ""}
       if [ -f $(echo '~{select_first([assemblyFilename, "generated.assembled.bam"])}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{select_first([assemblyFilename, "generated.assembled.bam"])}' | sed 's/\.[^.]*$//').bai $(echo '~{select_first([assemblyFilename, "generated.assembled.bam"])}' ).bai; fi
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 8, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "gridss/gridss:v2.2.3"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 31, 4])}G"
       preemptible: 2
     }
     output {
       File vcf = select_first([outputFilename, "generated.vcf"])
       File assembly = select_first([assemblyFilename, "generated.assembled.bam"])
       File assembly_bai = select_first([assemblyFilename, "generated.assembled.bam"]) + ".bai"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Gridss
   doc: |
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

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: gridss/gridss:v2.2.3

   inputs:
   - id: outputFilename
     label: outputFilename
     doc: (O=) VCF structural variation calls. Required.
     type:
     - string
     - 'null'
     default: generated.vcf
     inputBinding:
       prefix: OUTPUT=
       separate: false
   - id: reference
     label: reference
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
       prefix: REFERENCE_SEQUENCE=
       separate: false
   - id: bams
     label: bams
     doc: |-
       (I=File Coordinate-sorted input BAM file. Default value: null. This option may be specified 0 or more times.
     type:
       type: array
       inputBinding:
         prefix: INPUT=
         separate: false
       items: File
     inputBinding: {}
   - id: assemblyFilename
     label: assemblyFilename
     doc: Breakend assemblies which have undergone split read identification Required.
     type:
     - string
     - 'null'
     default: generated.assembled.bam
     inputBinding:
       prefix: ASSEMBLY=
       separate: false
   - id: inputLabel
     label: inputLabel
     doc: |-
       Input label. Variant calling evidence breakdowns are reported for each label. Default labels correspond to INPUT filenames. When specifying labels, labels must be provided for all input files. Default value: null. This option may be specified 0 or more times.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: INPUT_LABEL=
       separate: false
   - id: inputMaxFragmentSize
     label: inputMaxFragmentSize
     doc: |-
       Per input maximum concordant fragment size. Default value: null. This option may be specified 0 or more times.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: INPUT_MAX_FRAGMENT_SIZE=
       separate: false
   - id: inputMinFragmentSize
     label: inputMinFragmentSize
     doc: |-
       Per input minimum concordant fragment size. Default value: null. This option may be specified 0 or more times.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: INPUT_MIN_FRAGMENT_SIZE=
       separate: false
   - id: readPairConcordantPercent
     label: readPairConcordantPercent
     doc: |-
       Percent of read pairs considered concorant (0.0-1.0). If this is unset, the SAM proper pair flag is used to determine whether a read is discordantly aligned. Explicit fragment size specification overrides this setting. Default value: 0.995. This option can be set to 'null' to clear the default value.
     type:
     - float
     - 'null'
     inputBinding:
       prefix: READ_PAIR_CONCORDANT_PERCENT=
       separate: false
   - id: blacklist
     label: blacklist
     doc: |-
       (BL=File) BED blacklist of regions to ignore. Assembly of regions such as high-coverage centromeric repeats is slow, and if such regions are to be filtered in downstream analysis anyway, blacklisting those region will improve runtime performance. For human WGS, the ENCODE DAC blacklist is recommended. Default value: null.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: BLACKLIST=
       separate: false
   - id: configurationFile
     label: configurationFile
     doc: '(C=File) gridss configuration file containing overrides Default value: null.'
     type:
     - File
     - 'null'
     inputBinding:
       prefix: CONFIGURATION_FILE=
       separate: false
   - id: workerThreads
     label: workerThreads
     doc: |-
       (THREADS=Integer  Number of worker threads to spawn. Defaults to number of cores available. Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread. Default value: 6. This option can be set to 'null' to clear the default value.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: WORKER_THREADS=
       separate: false
   - id: workingDir
     label: workingDir
     doc: |-
       Directory to place intermediate results directories. Default location is the same directory as the associated input or output file. Default value: null.
     type: string
     default: .
     inputBinding:
       prefix: WORKING_DIR=
       separate: false
   - id: ignoreDuplicates
     label: ignoreDuplicates
     doc: |-
       Ignore reads marked as duplicates. Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: IGNORE_DUPLICATES=
       separate: false

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
     secondaryFiles:
     - |-
       ${

               function resolveSecondary(base, secPattern) {
                 if (secPattern[0] == "^") {
                   var spl = base.split(".");
                   var endIndex = spl.length > 1 ? spl.length - 1 : 1;
                   return resolveSecondary(spl.slice(undefined, endIndex).join("."), secPattern.slice(1));
                 }
                 return base + secPattern
               }
               return [
                       {
                           path: resolveSecondary(self.path, "^.bai"),
                           basename: resolveSecondary(self.basename, ".bai"),
                           class: "File",
                       }
               ];

       }
     outputBinding:
       glob: generated.assembled.bam
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - java
   - -XX:+UnlockExperimentalVMOptions
   - -XX:+UseCGroupMemoryLimitForHeap
   - -XX:MaxRAMFraction=1
   - -XshowSettings:vm
   - -Dsamjdk.create_index=true
   - -Dsamjdk.use_async_io_read_samtools=true
   - -Dsamjdk.use_async_io_write_samtools=true
   - -Dsamjdk.use_async_io_write_tribble=true
   - -Dgridss.gridss.output_to_temp_file=true
   - -Dsamjdk.buffer_size=4194304
   - -cp
   - /data/gridss/gridss-2.2.3-gridss-jar-with-dependencies.jar
   - gridss.CallVariants
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: gridss


