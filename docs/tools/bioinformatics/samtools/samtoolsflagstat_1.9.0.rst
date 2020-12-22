:orphan:

SamTools: Flagstat
=====================================

``SamToolsFlagstat`` · *1 contributor · 2 versions*

Does a full pass through the input file to calculate and print statistics to stdout.

Provides counts for each of 13 categories based primarily on bit flags in the FLAG field. Each category in the output is broken down into QC pass and QC fail. In the default output format, these are presented as "#PASS + #FAIL" followed by a description of the category.

The first row of output gives the total number of reads that are QC pass and fail (according to flag bit 0x200). For example:

122 + 28 in total (QC-passed reads + QC-failed reads)

Which would indicate that there are a total of 150 reads in the input file, 122 of which are marked as QC pass and 28 of which are marked as "not passing quality controls"

Following this, additional categories are given for reads which are:

secondary     0x100 bit set

supplementary     0x800 bit set

duplicates     0x400 bit set

mapped     0x4 bit not set

paired in sequencing     0x1 bit set

read1     both 0x1 and 0x40 bits set

read2     both 0x1 and 0x80 bits set

properly paired     both 0x1 and 0x2 bits set and 0x4 bit not set

with itself and mate mapped     0x1 bit set and neither 0x4 nor 0x8 bits set

singletons     both 0x1 and 0x8 bits set and bit 0x4 not set

And finally, two rows are given that additionally filter on the reference name (RNAME), mate reference name (MRNM), and mapping quality (MAPQ) fields:

with mate mapped to a different chr     0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME

with mate mapped to a different chr (mapQ>=5)     0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME and MAPQ >= 5)


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.samtools.flagstat.flagstat import SamToolsFlagstat_1_9

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "samtoolsflagstat_step",
           SamToolsFlagstat_1_9(
               bam=None,
           )
       )
       wf.output("out", source=samtoolsflagstat_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SamToolsFlagstat:

.. code-block:: bash

   # user inputs
   janis inputs SamToolsFlagstat > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam




5. Run SamToolsFlagstat with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SamToolsFlagstat





Information
------------

:ID: ``SamToolsFlagstat``
:URL: `http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS <http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS>`_
:Versions: 1.9.0, 1.7.0
:Container: quay.io/biocontainers/samtools:1.9--h8571acd_11
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-02-14
:Updated: 2020-02-14


Outputs
-----------

======  ================  ===============
name    type              documentation
======  ================  ===============
out     stdout<TextFile>
======  ================  ===============


Additional configuration (inputs)
---------------------------------

=======  =================  ========  ==========  ========================================================================
name     type               prefix      position  documentation
=======  =================  ========  ==========  ========================================================================
bam      BAM                                  10
threads  Optional<Integer>  -@                 5  Number of BAM compression threads to use in addition to main thread [0].
=======  =================  ========  ==========  ========================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task SamToolsFlagstat {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File bam
       Int? threads
     }
     command <<<
       set -e
       samtools flagstat \
         ~{if defined(threads) then ("-@ " + threads) else ''} \
         '~{bam}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/samtools:1.9--h8571acd_11"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = stdout()
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'SamTools: Flagstat'
   doc: |-
     Does a full pass through the input file to calculate and print statistics to stdout.

     Provides counts for each of 13 categories based primarily on bit flags in the FLAG field. Each category in the output is broken down into QC pass and QC fail. In the default output format, these are presented as "#PASS + #FAIL" followed by a description of the category.

     The first row of output gives the total number of reads that are QC pass and fail (according to flag bit 0x200). For example:

     122 + 28 in total (QC-passed reads + QC-failed reads)

     Which would indicate that there are a total of 150 reads in the input file, 122 of which are marked as QC pass and 28 of which are marked as "not passing quality controls"

     Following this, additional categories are given for reads which are:

     secondary     0x100 bit set

     supplementary     0x800 bit set

     duplicates     0x400 bit set

     mapped     0x4 bit not set

     paired in sequencing     0x1 bit set

     read1     both 0x1 and 0x40 bits set

     read2     both 0x1 and 0x80 bits set

     properly paired     both 0x1 and 0x2 bits set and 0x4 bit not set

     with itself and mate mapped     0x1 bit set and neither 0x4 nor 0x8 bits set

     singletons     both 0x1 and 0x8 bits set and bit 0x4 not set

     And finally, two rows are given that additionally filter on the reference name (RNAME), mate reference name (MRNM), and mapping quality (MAPQ) fields:

     with mate mapped to a different chr     0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME

     with mate mapped to a different chr (mapQ>=5)     0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME and MAPQ >= 5)

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/samtools:1.9--h8571acd_11

   inputs:
   - id: bam
     label: bam
     type: File
     inputBinding:
       position: 10
   - id: threads
     label: threads
     doc: Number of BAM compression threads to use in addition to main thread [0].
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -@
       position: 5

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - samtools
   - flagstat
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: SamToolsFlagstat


