:orphan:

Bcl2Fastq
=====================

``bcl2fastq`` · *1 contributor · 1 version*

BCL to FASTQ file converter

.. warning::

   Bcl2Fastq did not include a container. You can provide one through the command line by including
   the following instruction:

   .. code-block:: bash

      janis run --container-override 'bcl2fastq=<organisation/container:version>' bcl2fastq
    
Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.illumina.bcl2fastq.versions import Bcl2Fastq_2_20_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bcl2fastq_step",
           Bcl2Fastq_2_20_0(
               runFolderDir=None,
               sampleSheet=None,
               loadingThreads=None,
               processingThreads=None,
               writingThreads=None,
           )
       )
       wf.output("unalignedReads", source=bcl2fastq_step.unalignedReads)
       wf.output("stats", source=bcl2fastq_step.stats)
       wf.output("interop", source=bcl2fastq_step.interop)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bcl2fastq:

.. code-block:: bash

   # user inputs
   janis inputs bcl2fastq > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       runFolderDir: null
       sampleSheet: sampleSheet.csv




5. Run bcl2fastq with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       --container-override 'bcl2fastq=<organisation/container:version>' \
       bcl2fastq





Information
------------

:ID: ``bcl2fastq``
:URL: `https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html <https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html>`_
:Versions: 2.20.0
:Container: None
:Authors: Matthias De Smet (@mattdsm)
:Citations: None
:Created: 2020-03-05
:Updated: 2020-03-05


Outputs
-----------

==============  ==============  ===============
name            type            documentation
==============  ==============  ===============
unalignedReads  Array<FastqGz>
stats           Array<File>
interop         Array<File>
==============  ==============  ===============


Additional configuration (inputs)
---------------------------------

=============================  =================  ===================================  ==========  ===============================================================================================================================
name                           type               prefix                               position    documentation
=============================  =================  ===================================  ==========  ===============================================================================================================================
runFolderDir                   Directory          -R                                               path to runfolder directory
sampleSheet                    csv                --sample-sheet                                   path to the sample sheet
loadingThreads                 Integer            -r                                               number of threads used for loading BCL data
processingThreads              Integer            -p                                               number of threads used for processing demultiplexed data
writingThreads                 Integer            -w                                               number of threads used for writing FASTQ data
minimumTrimmedReadLength       Optional<Integer>  --minimum-trimmed-read-length                    minimum read length after adapter trimming
useBasesMask                   Optional<String>   --use-bases-mask                                 specifies how to use each cycle
maskShortAdapterReads          Optional<Integer>  --mask-short-adapter-reads                       smallest number of remaining bases (after masking bases below the minimum trimmed read length) below which whole read is masked
adapterStringency              Optional<Float>    --adapter-stringency                             adapter stringency
ignoreMissingBcls              Optional<Boolean>  --ignore-missing-bcls                            assume 'N'/'#' for missing calls
ignoreMissingFilter            Optional<Boolean>  --ignore-missing-filter                          assume 'true' for missing filters
ignoreMissingPositions         Optional<Boolean>  --ignore-missing-positions                       assume [0,i] for missing positions, where i is incremented starting from 0
writeFastqReverseComplement    Optional<Boolean>  --write-fastq-reverse-complement                 generate FASTQs containing reverse complements of actual data
withFailedReads                Optional<Boolean>  --with-failed-reads                              include non-PF clusters
createFastqForIndexReads       Optional<Boolean>  --create-fastq-for-index-reads                   create FASTQ files also for index reads
findAdaptersWithSlidingWindow  Optional<Boolean>  --find-adapters-with-sliding-window              find adapters with simple sliding window algorithm
noBgzfCompression              Optional<Boolean>  --no-bgzf-compression                            turn off BGZF compression for FASTQ files
barcodeMismatches              Optional<Integer>  --barcode-mismatches                             number of allowed mismatches per index
noLaneSplitting                Optional<Boolean>  --no-lane-splitting                              do not split fastq files by lane
=============================  =================  ===================================  ==========  ===============================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task bcl2fastq {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Directory runFolderDir
       File sampleSheet
       Int? loadingThreads
       Int? processingThreads
       Int? writingThreads
       Int? minimumTrimmedReadLength
       String? useBasesMask
       Int? maskShortAdapterReads
       Float? adapterStringency
       Boolean? ignoreMissingBcls
       Boolean? ignoreMissingFilter
       Boolean? ignoreMissingPositions
       Boolean? writeFastqReverseComplement
       Boolean? withFailedReads
       Boolean? createFastqForIndexReads
       Boolean? findAdaptersWithSlidingWindow
       Boolean? noBgzfCompression
       Int? barcodeMismatches
       Boolean? noLaneSplitting
     }
     command <<<
       set -e
       bcl2fastq \
         -R '~{runFolderDir}' \
         --sample-sheet '~{sampleSheet}' \
         -r ~{select_first([loadingThreads, 4])} \
         -p ~{select_first([processingThreads, 4])} \
         -w ~{select_first([writingThreads, 4])} \
         ~{if defined(minimumTrimmedReadLength) then ("--minimum-trimmed-read-length " + minimumTrimmedReadLength) else ''} \
         ~{if defined(useBasesMask) then ("--use-bases-mask '" + useBasesMask + "'") else ""} \
         ~{if defined(maskShortAdapterReads) then ("--mask-short-adapter-reads " + maskShortAdapterReads) else ''} \
         ~{if defined(adapterStringency) then ("--adapter-stringency " + adapterStringency) else ''} \
         ~{if (defined(ignoreMissingBcls) && select_first([ignoreMissingBcls])) then "--ignore-missing-bcls" else ""} \
         ~{if (defined(ignoreMissingFilter) && select_first([ignoreMissingFilter])) then "--ignore-missing-filter" else ""} \
         ~{if (defined(ignoreMissingPositions) && select_first([ignoreMissingPositions])) then "--ignore-missing-positions" else ""} \
         ~{if (defined(writeFastqReverseComplement) && select_first([writeFastqReverseComplement])) then "--write-fastq-reverse-complement" else ""} \
         ~{if (defined(withFailedReads) && select_first([withFailedReads])) then "--with-failed-reads" else ""} \
         ~{if (defined(createFastqForIndexReads) && select_first([createFastqForIndexReads])) then "--create-fastq-for-index-reads" else ""} \
         ~{if (defined(findAdaptersWithSlidingWindow) && select_first([findAdaptersWithSlidingWindow])) then "--find-adapters-with-sliding-window" else ""} \
         ~{if (defined(noBgzfCompression) && select_first([noBgzfCompression])) then "--no-bgzf-compression" else ""} \
         ~{if defined(barcodeMismatches) then ("--barcode-mismatches " + barcodeMismatches) else ''} \
         ~{if (defined(noLaneSplitting) && select_first([noLaneSplitting])) then " --no-lane-splitting" else ""} \
         --output-dir '.'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4, 4])}G"
       preemptible: 2
     }
     output {
       Array[File] unalignedReads = glob("*/*.fastq.gz")
       Array[File] stats = glob("Stats/*")
       Array[File] interop = glob("InterOp/*")
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Bcl2Fastq
   doc: BCL to FASTQ file converter

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement

   inputs:
   - id: runFolderDir
     label: runFolderDir
     doc: path to runfolder directory
     type: Directory
     inputBinding:
       prefix: -R
   - id: sampleSheet
     label: sampleSheet
     doc: path to the sample sheet
     type: File
     inputBinding:
       prefix: --sample-sheet
   - id: loadingThreads
     label: loadingThreads
     doc: number of threads used for loading BCL data
     type: int
     default: 4
     inputBinding:
       prefix: -r
   - id: processingThreads
     label: processingThreads
     doc: number of threads used for processing demultiplexed data
     type: int
     default: 4
     inputBinding:
       prefix: -p
   - id: writingThreads
     label: writingThreads
     doc: number of threads used for writing FASTQ data
     type: int
     default: 4
     inputBinding:
       prefix: -w
   - id: minimumTrimmedReadLength
     label: minimumTrimmedReadLength
     doc: minimum read length after adapter trimming
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --minimum-trimmed-read-length
   - id: useBasesMask
     label: useBasesMask
     doc: specifies how to use each cycle
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --use-bases-mask
   - id: maskShortAdapterReads
     label: maskShortAdapterReads
     doc: |-
       smallest number of remaining bases (after masking bases below the minimum trimmed read length) below which whole read is masked
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --mask-short-adapter-reads
   - id: adapterStringency
     label: adapterStringency
     doc: adapter stringency
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --adapter-stringency
   - id: ignoreMissingBcls
     label: ignoreMissingBcls
     doc: assume 'N'/'#' for missing calls
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ignore-missing-bcls
   - id: ignoreMissingFilter
     label: ignoreMissingFilter
     doc: assume 'true' for missing filters
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ignore-missing-filter
   - id: ignoreMissingPositions
     label: ignoreMissingPositions
     doc: assume [0,i] for missing positions, where i is incremented starting from 0
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ignore-missing-positions
   - id: writeFastqReverseComplement
     label: writeFastqReverseComplement
     doc: generate FASTQs containing reverse complements of actual data
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --write-fastq-reverse-complement
   - id: withFailedReads
     label: withFailedReads
     doc: include non-PF clusters
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --with-failed-reads
   - id: createFastqForIndexReads
     label: createFastqForIndexReads
     doc: create FASTQ files also for index reads
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --create-fastq-for-index-reads
   - id: findAdaptersWithSlidingWindow
     label: findAdaptersWithSlidingWindow
     doc: find adapters with simple sliding window algorithm
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --find-adapters-with-sliding-window
   - id: noBgzfCompression
     label: noBgzfCompression
     doc: turn off BGZF compression for FASTQ files
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-bgzf-compression
   - id: barcodeMismatches
     label: barcodeMismatches
     doc: number of allowed mismatches per index
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --barcode-mismatches
   - id: noLaneSplitting
     label: noLaneSplitting
     doc: do not split fastq files by lane
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: ' --no-lane-splitting'

   outputs:
   - id: unalignedReads
     label: unalignedReads
     type:
       type: array
       items: File
     outputBinding:
       glob: '*/*.fastq.gz'
       loadContents: false
   - id: stats
     label: stats
     type:
       type: array
       items: File
     outputBinding:
       glob: Stats/*
       loadContents: false
   - id: interop
     label: interop
     type:
       type: array
       items: File
     outputBinding:
       glob: InterOp/*
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: bcl2fastq
   arguments:
   - prefix: --output-dir
     position: 0
     valueFrom: .

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: bcl2fastq


