:orphan:

Trimmomatic: Paired End (PE)
===================================================

``trimmomaticPairedEnd`` · *1 contributor · 1 version*

Trimmomatic is a fast, multithreaded command line tool that can be used to trim and crop
Illumina (FASTQ) data as well as to remove adapters. These adapters can pose a real problem
depending on the library preparation and downstream application.

There are two major modes of the program: Paired end mode and Single end mode. The
paired end mode will maintain correspondence of read pairs and also use the additional
information contained in paired reads to better find adapter or PCR primer fragments
introduced by the library preparation process.

Trimmomatic works with FASTQ files (using phred + 33 or phred + 64 quality scores,
depending on the Illumina pipeline used). Files compressed using either "gzip" or "bzip2" are
supported, and are identified by use of ".gz" or ".bz2" file extensions. 


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.usadellab.trimmomatic.versions import TrimmomaticPairedEnd_0_35

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "trimmomaticpairedend_step",
           TrimmomaticPairedEnd_0_35(
               steps=None,
               sampleName=None,
               inp=None,
           )
       )
       wf.output("pairedOut", source=trimmomaticpairedend_step.pairedOut)
       wf.output("unpairedOut", source=trimmomaticpairedend_step.unpairedOut)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for trimmomaticPairedEnd:

.. code-block:: bash

   # user inputs
   janis inputs trimmomaticPairedEnd > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inp:
       - inp_0.fastq.gz
       - inp_1.fastq.gz
       sampleName: <value>
       steps:
       - <value>
       - <value>




5. Run trimmomaticPairedEnd with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       trimmomaticPairedEnd





Information
------------

:ID: ``trimmomaticPairedEnd``
:URL: `http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf <http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf>`_
:Versions: 0.35
:Container: quay.io/biocontainers/trimmomatic:0.35--6
:Authors: illusional
:Citations: Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.
:DOI: 10.1093/bioinformatics/btu170
:Created: 2020-05-25
:Updated: 2020-05-25


Outputs
-----------

===========  ===========  ===============
name         type         documentation
===========  ===========  ===============
pairedOut    FastqGzPair
unpairedOut  FastqGzPair
===========  ===========  ===============


Additional configuration (inputs)
---------------------------------

=========================  ==================  ========  ==========  ======================================================================================================
name                       type                prefix      position  documentation
=========================  ==================  ========  ==========  ======================================================================================================
steps                      Array<String>                        100  ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
                                                                     SLIDINGWINDOW: Performs a sliding window trimming approach. It starts
                                                                     scanning at the 5" end and clips the read once the average quality within the window
                                                                     falls below a threshold.
                                                                     MAXINFO: An adaptive quality trimmer which balances read length and error rate to
                                                                     maximise the value of each read
                                                                     LEADING: Cut bases off the start of a read, if below a threshold quality
                                                                     TRAILING: Cut bases off the end of a read, if below a threshold quality
                                                                     CROP: Cut the read to a specified length by removing bases from the end
                                                                     HEADCROP: Cut the specified number of bases from the start of the read
                                                                     MINLEN: Drop the read if it is below a specified length
                                                                     AVGQUAL: Drop the read if the average quality is below the specified level
                                                                     TOPHRED33: Convert quality scores to Phred-33
                                                                     TOPHRED64: Convert quality scores to Phred-64
sampleName                 String                                    Used to name the output
inp                        FastqGzPair                            5
threads                    Optional<Integer>   -threads           2
phred33                    Optional<Boolean>   -phred33           3  Use phred + 33 quality score. If no quality encoding is specified, it will be determined automatically
phred64                    Optional<Boolean>   -phred64           3  Use phred + 64 quality score. If no quality encoding is specified, it will be determined automatically
trimLogFilename            Optional<Filename>  -trimlog           4  Specifying a trimlog file creates a log of all read trimmings, indicating the following details:

                                                                         - the read name
                                                                         - the surviving sequence length
                                                                         - the location of the first surviving base, aka. the amount trimmed from the start
                                                                         - the location of the last surviving base in the original read
                                                                         - the amount trimmed from the end
outputFilename_R1          Optional<Filename>                     6
outputFilenameUnpaired_R1  Optional<Filename>                     7
outputFilename_R2          Optional<Filename>                     8
outputFilenameUnpaired_R2  Optional<Filename>                     9
=========================  ==================  ========  ==========  ======================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task trimmomaticPairedEnd {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String] steps
       String sampleName
       Int? threads
       Boolean? phred33
       Boolean? phred64
       String? trimLogFilename
       Array[File] inp
       String? outputFilename_R1
       String? outputFilenameUnpaired_R1
       String? outputFilename_R2
       String? outputFilenameUnpaired_R2
     }
     command <<<
       set -e
       trimmomatic \
         'PE' \
         ~{if defined(threads) then ("-threads " + threads) else ''} \
         ~{if (defined(phred33) && select_first([phred33])) then "-phred33" else ""} \
         ~{if (defined(phred64) && select_first([phred64])) then "-phred64" else ""} \
         -trimlog '~{select_first([trimLogFilename, "trimlog.log"])}' \
         ~{if length(inp) > 0 then "'" + sep("' '", inp) + "'" else ""} \
         '~{select_first([outputFilename_R1, "~{sampleName}-_R1.trimmed.fastq.gz"])}' \
         '~{select_first([outputFilenameUnpaired_R1, "~{sampleName}-_R1.unpaired.fastq.gz"])}' \
         '~{select_first([outputFilename_R2, "~{sampleName}-_R2.trimmed.fastq.gz"])}' \
         '~{select_first([outputFilenameUnpaired_R2, "~{sampleName}-_R2.unpaired.fastq.gz"])}' \
         ~{if length(steps) > 0 then "'" + sep("' '", steps) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/trimmomatic:0.35--6"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       Array[File] pairedOut = glob("*trimmed.fastq.gz")
       Array[File] unpairedOut = glob("*unpaired.fastq.gz")
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'Trimmomatic: Paired End (PE)'
   doc: |-
     Trimmomatic is a fast, multithreaded command line tool that can be used to trim and crop
     Illumina (FASTQ) data as well as to remove adapters. These adapters can pose a real problem
     depending on the library preparation and downstream application.

     There are two major modes of the program: Paired end mode and Single end mode. The
     paired end mode will maintain correspondence of read pairs and also use the additional
     information contained in paired reads to better find adapter or PCR primer fragments
     introduced by the library preparation process.

     Trimmomatic works with FASTQ files (using phred + 33 or phred + 64 quality scores,
     depending on the Illumina pipeline used). Files compressed using either "gzip" or "bzip2" are
     supported, and are identified by use of ".gz" or ".bz2" file extensions. 

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/trimmomatic:0.35--6

   inputs:
   - id: steps
     label: steps
     doc: |
       ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
       SLIDINGWINDOW: Performs a sliding window trimming approach. It starts
       scanning at the 5" end and clips the read once the average quality within the window
       falls below a threshold.
       MAXINFO: An adaptive quality trimmer which balances read length and error rate to
       maximise the value of each read
       LEADING: Cut bases off the start of a read, if below a threshold quality
       TRAILING: Cut bases off the end of a read, if below a threshold quality
       CROP: Cut the read to a specified length by removing bases from the end
       HEADCROP: Cut the specified number of bases from the start of the read
       MINLEN: Drop the read if it is below a specified length
       AVGQUAL: Drop the read if the average quality is below the specified level
       TOPHRED33: Convert quality scores to Phred-33
       TOPHRED64: Convert quality scores to Phred-64
     type:
       type: array
       items: string
     inputBinding:
       position: 100
   - id: sampleName
     label: sampleName
     doc: Used to name the output
     type: string
   - id: threads
     label: threads
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -threads
       position: 2
   - id: phred33
     label: phred33
     doc: |-
       Use phred + 33 quality score. If no quality encoding is specified, it will be determined automatically
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -phred33
       position: 3
   - id: phred64
     label: phred64
     doc: |-
       Use phred + 64 quality score. If no quality encoding is specified, it will be determined automatically
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -phred64
       position: 3
   - id: trimLogFilename
     label: trimLogFilename
     doc: |-
       Specifying a trimlog file creates a log of all read trimmings, indicating the following details:

           - the read name
           - the surviving sequence length
           - the location of the first surviving base, aka. the amount trimmed from the start
           - the location of the last surviving base in the original read
           - the amount trimmed from the end
     type:
     - string
     - 'null'
     default: trimlog.log
     inputBinding:
       prefix: -trimlog
       position: 4
   - id: inp
     label: inp
     type:
       type: array
       items: File
     inputBinding:
       position: 5
       itemSeparator: ' '
   - id: outputFilename_R1
     label: outputFilename_R1
     type:
     - string
     - 'null'
     default: generated-_R1.trimmed.fastq.gz
     inputBinding:
       position: 6
       valueFrom: $(inputs.sampleName)-_R1.trimmed.fastq.gz
   - id: outputFilenameUnpaired_R1
     label: outputFilenameUnpaired_R1
     type:
     - string
     - 'null'
     default: generated-_R1.unpaired.fastq.gz
     inputBinding:
       position: 7
       valueFrom: $(inputs.sampleName)-_R1.unpaired.fastq.gz
   - id: outputFilename_R2
     label: outputFilename_R2
     type:
     - string
     - 'null'
     default: generated-_R2.trimmed.fastq.gz
     inputBinding:
       position: 8
       valueFrom: $(inputs.sampleName)-_R2.trimmed.fastq.gz
   - id: outputFilenameUnpaired_R2
     label: outputFilenameUnpaired_R2
     type:
     - string
     - 'null'
     default: generated-_R2.unpaired.fastq.gz
     inputBinding:
       position: 9
       valueFrom: $(inputs.sampleName)-_R2.unpaired.fastq.gz

   outputs:
   - id: pairedOut
     label: pairedOut
     type:
       type: array
       items: File
     outputBinding:
       glob: '*trimmed.fastq.gz'
       loadContents: false
   - id: unpairedOut
     label: unpairedOut
     type:
       type: array
       items: File
     outputBinding:
       glob: '*unpaired.fastq.gz'
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: trimmomatic
   arguments:
   - position: 0
     valueFrom: PE

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: trimmomaticPairedEnd


