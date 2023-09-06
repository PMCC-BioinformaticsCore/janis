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

       from janis_bioinformatics.tools.samtools.flagstat.flagstat import SamToolsFlagstat_1_7

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "samtoolsflagstat_step",
           SamToolsFlagstat_1_7(
               bam=None,
           )
       )
       wf.output("out", source=samtoolsflagstat_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

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

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          SamToolsFlagstat

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``SamToolsFlagstat``
:URL: `http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS <http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS>`_
:Versions: 1.9.0, 1.7.0
:Container: biocontainers/samtools:v1.7.0_cv3
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-02-14
:Updated: 2020-02-14


Outputs
-----------

======  ========  ===============
name    type      documentation
======  ========  ===============
out     TextFile
======  ========  ===============


Additional configuration (inputs)
---------------------------------

==============  ==================  ========  ==========  ========================================================================
name            type                prefix      position  documentation
==============  ==================  ========  ==========  ========================================================================
bam             BAM                                   10
threads         Optional<Integer>   -@                 5  Number of BAM compression threads to use in addition to main thread [0].
outputFilename  Optional<Filename>  >                 11
==============  ==================  ========  ==========  ========================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task SamToolsFlagstat {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       File bam
       Int? threads
       String? outputFilename
     }

     command <<<
       set -e
       samtools flagstat \
         ~{if defined(threads) then ("-@ " + threads) else ''} \
         '~{bam}' \
         > '~{select_first([outputFilename, "generated"])}'
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "biocontainers/samtools:v1.7.0_cv3"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }

     output {
       File out = select_first([outputFilename, "generated"])
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'SamTools: Flagstat'

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: biocontainers/samtools:v1.7.0_cv3

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
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: '>'
       position: 11

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated
       loadContents: false
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


