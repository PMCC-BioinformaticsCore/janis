:orphan:

STAR Aligner
===========================

*1 contributor · 1 version*

Usage: STAR  [options]... --genomeDir /path/to/genome/index/   --readFilesIn R1.fq R2.fq
Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2019

For more details see:
<https://github.com/alexdobin/STAR>
<https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>
            


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.star.versions import StarAligner_2_7_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "star_aligner_step",
           StarAligner_2_7_1(

           )
       )
       wf.output("logFinalOut", source=star_aligner_step.logFinalOut)
       wf.output("logOut", source=star_aligner_step.logOut)
       wf.output("logProgressOut", source=star_aligner_step.logProgressOut)
       wf.output("sjOutTab", source=star_aligner_step.sjOutTab)
       wf.output("out", source=star_aligner_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for star_aligner:

.. code-block:: bash

   # user inputs
   janis inputs star_aligner > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run star_aligner with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       star_aligner





Information
------------

:ID: ``star_aligner``
:URL: `https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>`_
:Versions: v2.7.1a
:Container: quay.io/biocontainers/star:2.7.3a--0
:Authors: Jiaan Yu
:Citations: Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics (Oxford, England), 29(1), 15–21. https://doi.org/10.1093/bioinformatics/bts635
:DOI: https://doi.org/10.1093/bioinformatics/bts635
:Created: 2020-04-02
:Updated: 2020-04-02


Outputs
-----------

==============  ======  ===============
name            type    documentation
==============  ======  ===============
logFinalOut     File
logOut          File
logProgressOut  File
sjOutTab        File
out             BAM
==============  ======  ===============


Additional configuration (inputs)
---------------------------------

=================  ========================  ===================  ==========  ========================================================================================================================================================================================================================================================================================================================
name               type                      prefix               position    documentation
=================  ========================  ===================  ==========  ========================================================================================================================================================================================================================================================================================================================
help               Optional<Boolean>         --help                           help page
runThreadN         Optional<Integer>         --runThreadN                     int: number of threads to run STAR. Default: 1.
genomeDir          Optional<Directory>       --genomeDir                      string: path to the directory where genome files are stored (for –runMode alignReads) or will be generated (for –runMode generateGenome). Default: ./GenomeDir
readFilesIn        Optional<Array<FastqGz>>  --readFilesIn                    string(s): paths to files that contain input read1 (and, if needed, read2). Default: Read1,Read2.
outFileNamePrefix  Optional<Filename>        --outFileNamePrefix              string: output files name prefix (including full or relative path). Can only be defined on the command line.
outSAMtype         Optional<Array<String>>   --outSAMtype                     strings: type of SAM/BAM output. 1st word: "BAM": outputBAMwithoutsorting, "SAM": outputSAMwithoutsorting, "None": no SAM/BAM output. 2nd,3rd: "Unsorted": standard unsorted. "SortedByCoordinate": sorted by coordinate. This option will allocate extra memory for sorting which can be specified by –limitBAMsortRAM.
outSAMunmapped     Optional<String>          --outSAMunmapped                 string(s): output of unmapped reads in the SAM format
outSAMattributes   Optional<String>          --outSAMattributes               string: a string of desired SAM attributes, in the order desired for the output SAM
readFilesCommand   Optional<String>          --readFilesCommand               string(s): command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout
=================  ========================  ===================  ==========  ========================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task star_aligner {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Boolean? help
       Int? runThreadN
       Directory? genomeDir
       Array[File]? readFilesIn
       String? outFileNamePrefix
       Array[String]? outSAMtype
       String? outSAMunmapped
       String? outSAMattributes
       String? readFilesCommand
     }
     command <<<
       set -e
       STAR \
         ~{if defined(help) then "--help" else ""} \
         ~{if defined(select_first([runThreadN, select_first([runtime_cpu, 1])])) then ("--runThreadN " + select_first([runThreadN, select_first([runtime_cpu, 1])])) else ''} \
         ~{if defined(genomeDir) then ("--genomeDir '" + genomeDir + "'") else ""} \
         ~{if (defined(readFilesIn) && length(select_first([readFilesIn])) > 0) then "--readFilesIn '" + sep("','", select_first([readFilesIn])) + "'" else ""} \
         --outFileNamePrefix '~{select_first([outFileNamePrefix, "generated"])}' \
         ~{if (defined(outSAMtype) && length(select_first([outSAMtype])) > 0) then "--outSAMtype '" + sep("' '", select_first([outSAMtype])) + "'" else ""} \
         ~{if defined(outSAMunmapped) then ("--outSAMunmapped '" + outSAMunmapped + "'") else ""} \
         ~{if defined(outSAMattributes) then ("--outSAMattributes '" + outSAMattributes + "'") else ""} \
         ~{if defined(readFilesCommand) then ("--readFilesCommand '" + readFilesCommand + "'") else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/star:2.7.3a--0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 32, 4])}G"
       preemptible: 2
     }
     output {
       File logFinalOut = (select_first([outFileNamePrefix, "generated"]) + "Log.final.out")
       File logOut = (select_first([outFileNamePrefix, "generated"]) + "Log.out")
       File logProgressOut = (select_first([outFileNamePrefix, "generated"]) + "Log.progress.out")
       File sjOutTab = (select_first([outFileNamePrefix, "generated"]) + "SJ.out.tab")
       File out = (select_first([outFileNamePrefix, "generated"]) + "Aligned.out.bam")
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.0
   label: STAR Aligner
   doc: |-
     Usage: STAR  [options]... --genomeDir /path/to/genome/index/   --readFilesIn R1.fq R2.fq
     Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2019

     For more details see:
     <https://github.com/alexdobin/STAR>
     <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>
              

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/star:2.7.3a--0

   inputs:
   - id: help
     label: help
     doc: help page
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --help
   - id: runThreadN
     label: runThreadN
     doc: 'int: number of threads to run STAR. Default: 1.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --runThreadN
       valueFrom: |-
         $([inputs.runtime_cpu, 4, 1].filter(function (inner) { return inner != null })[0])
   - id: genomeDir
     label: genomeDir
     doc: |-
       string: path to the directory where genome files are stored (for –runMode alignReads) or will be generated (for –runMode generateGenome). Default: ./GenomeDir
     type:
     - Directory
     - 'null'
     inputBinding:
       prefix: --genomeDir
   - id: readFilesIn
     label: readFilesIn
     doc: |-
       string(s): paths to files that contain input read1 (and, if needed, read2). Default: Read1,Read2.
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       prefix: --readFilesIn
       itemSeparator: ','
   - id: outFileNamePrefix
     label: outFileNamePrefix
     doc: |-
       string: output files name prefix (including full or relative path). Can only be defined on the command line.
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --outFileNamePrefix
   - id: outSAMtype
     label: outSAMtype
     doc: |-
       strings: type of SAM/BAM output. 1st word: "BAM": outputBAMwithoutsorting, "SAM": outputSAMwithoutsorting, "None": no SAM/BAM output. 2nd,3rd: "Unsorted": standard unsorted. "SortedByCoordinate": sorted by coordinate. This option will allocate extra memory for sorting which can be specified by –limitBAMsortRAM.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --outSAMtype
       itemSeparator: ' '
   - id: outSAMunmapped
     label: outSAMunmapped
     doc: 'string(s): output of unmapped reads in the SAM format'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outSAMunmapped
   - id: outSAMattributes
     label: outSAMattributes
     doc: |-
       string: a string of desired SAM attributes, in the order desired for the output SAM
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outSAMattributes
   - id: readFilesCommand
     label: readFilesCommand
     doc: |-
       string(s): command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --readFilesCommand

   outputs:
   - id: logFinalOut
     label: logFinalOut
     type: File
     outputBinding:
       glob: $((inputs.outFileNamePrefix + "Log.final.out"))
       outputEval: $((inputs.outFileNamePrefix + "Log.final.out"))
       loadContents: false
   - id: logOut
     label: logOut
     type: File
     outputBinding:
       glob: $((inputs.outFileNamePrefix + "Log.out"))
       outputEval: $((inputs.outFileNamePrefix + "Log.out"))
       loadContents: false
   - id: logProgressOut
     label: logProgressOut
     type: File
     outputBinding:
       glob: $((inputs.outFileNamePrefix + "Log.progress.out"))
       outputEval: $((inputs.outFileNamePrefix + "Log.progress.out"))
       loadContents: false
   - id: sjOutTab
     label: sjOutTab
     type: File
     outputBinding:
       glob: $((inputs.outFileNamePrefix + "SJ.out.tab"))
       outputEval: $((inputs.outFileNamePrefix + "SJ.out.tab"))
       loadContents: false
   - id: out
     label: out
     type: File
     outputBinding:
       glob: $((inputs.outFileNamePrefix + "Aligned.out.bam"))
       outputEval: $((inputs.outFileNamePrefix + "Aligned.out.bam"))
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: STAR
   arguments: []
   id: star_aligner


