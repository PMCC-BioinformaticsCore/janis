:orphan:

Bcl2Fastq
=====================

*1 contributor · 1 version*

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
       sampleSheet: sampleSheet




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
