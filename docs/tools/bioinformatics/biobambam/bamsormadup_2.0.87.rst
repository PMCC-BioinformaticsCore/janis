:orphan:

BamSorMaDup
=========================

``bamsormadup`` · *1 contributor · 1 version*

bamsormadup: parallel sorting and duplicate marking


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.biobambam.bamsormadup.versions import BamSorMaDup_2_0_87

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bamsormadup_step",
           BamSorMaDup_2_0_87(
               alignedReads=None,
           )
       )
       wf.output("out", source=bamsormadup_step.out)
       wf.output("metrics", source=bamsormadup_step.metrics)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bamsormadup:

.. code-block:: bash

   # user inputs
   janis inputs bamsormadup > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       alignedReads: alignedReads.bam




5. Run bamsormadup with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bamsormadup





Information
------------

:ID: ``bamsormadup``
:URL: `https://gitlab.com/german.tischler/biobambam2 <https://gitlab.com/german.tischler/biobambam2>`_
:Versions: 2.0.87
:Container: quay.io/biocontainers/biobambam:2.0.87--1
:Authors: Matthias De Smet (@mattdsm)
:Citations: None
:Created: 2020-02-26
:Updated: 2020-02-26


Outputs
-----------

=======  ===========  ===============
name     type         documentation
=======  ===========  ===============
out      stdout<BAM>
metrics  File
=======  ===========  ===============


Additional configuration (inputs)
---------------------------------

==============  ==================  ===============  ==========  =========================================================================================================
name            type                prefix             position  documentation
==============  ==================  ===============  ==========  =========================================================================================================
alignedReads    BAM                                         200
outputFilename  Optional<Filename>
level           Optional<Integer>   level=                       compression settings for output bam file (-1=zlib default,0=uncompressed,1=fast,9=best)
tempLevel       Optional<Integer>   templevel=                   compression settings for temporary bam files (-1=zlib default,0=uncompressed,1=fast,9=best)
threads         Optional<Integer>   threads=                     Number of threads. (default = 1)
sortOrder       Optional<String>    SO=                          output sort order(coordinate by default)
optMinPixelDif  Optional<Integer>   optminpixeldif=              pixel difference threshold for optical duplicates (patterned flowcell: 12000, unpatterned flowcell: 2500)
==============  ==================  ===============  ==========  =========================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task bamsormadup {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File alignedReads
       String? outputFilename
       Int? level
       Int? tempLevel
       Int? threads
       String? sortOrder
       Int? optMinPixelDif
     }
     command <<<
       set -e
       bamsormadup \
         ~{if defined(select_first([level, 0])) then ("level=" + select_first([level, 0])) else ''} \
         ~{if defined(select_first([tempLevel, 0])) then ("templevel=" + select_first([tempLevel, 0])) else ''} \
         ~{if defined(select_first([threads, select_first([runtime_cpu, 1])])) then ("threads=" + select_first([threads, select_first([runtime_cpu, 1])])) else ''} \
         ~{if defined(select_first([sortOrder, "coordinate"])) then ("SO='" + select_first([sortOrder, "coordinate"]) + "'") else ""} \
         ~{if defined(select_first([optMinPixelDif, 2500])) then ("optminpixeldif=" + select_first([optMinPixelDif, 2500])) else ''} \
         M= 'metrics.txt' \
         inputformat= 'bam' \
         outputFormat= 'bam' \
         '~{alignedReads}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/biobambam:2.0.87--1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 16, 4])}G"
       preemptible: 2
     }
     output {
       File out = stdout()
       File metrics = glob("metrics.txt")[0]
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: BamSorMaDup
   doc: 'bamsormadup: parallel sorting and duplicate marking'

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/biobambam:2.0.87--1

   inputs:
   - id: alignedReads
     label: alignedReads
     type: File
     inputBinding:
       position: 200
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.bam
   - id: level
     label: level
     doc: |-
       compression settings for output bam file (-1=zlib default,0=uncompressed,1=fast,9=best)
     type: int
     default: 0
     inputBinding:
       prefix: level=
       separate: false
   - id: tempLevel
     label: tempLevel
     doc: |-
       compression settings for temporary bam files (-1=zlib default,0=uncompressed,1=fast,9=best)
     type: int
     default: 0
     inputBinding:
       prefix: templevel=
       separate: false
   - id: threads
     label: threads
     doc: Number of threads. (default = 1)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: threads=
       valueFrom: |-
         $([inputs.runtime_cpu, 4, 1].filter(function (inner) { return inner != null })[0])
       separate: false
   - id: sortOrder
     label: sortOrder
     doc: output sort order(coordinate by default)
     type: string
     default: coordinate
     inputBinding:
       prefix: SO=
       separate: false
   - id: optMinPixelDif
     label: optMinPixelDif
     doc: |-
       pixel difference threshold for optical duplicates (patterned flowcell: 12000, unpatterned flowcell: 2500)
     type: int
     default: 2500
     inputBinding:
       prefix: optminpixeldif=
       separate: false

   outputs:
   - id: out
     label: out
     type: stdout
   - id: metrics
     label: metrics
     type: File
     outputBinding:
       glob: metrics.txt
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - bamsormadup
   arguments:
   - prefix: M=
     position: 0
     valueFrom: metrics.txt
     separate: false
   - prefix: inputformat=
     position: 0
     valueFrom: bam
     separate: false
   - prefix: outputFormat=
     position: 0
     valueFrom: bam
     separate: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: bamsormadup


