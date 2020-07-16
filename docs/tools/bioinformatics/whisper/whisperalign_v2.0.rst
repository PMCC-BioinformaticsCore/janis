:orphan:

Whisper-Align
============================

*1 contributor · 1 version*

Builds a whisper index


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.whisper.align.versions import WhisperAlign_2_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "whisperalign_step",
           WhisperAlign_2_0(
               index=None,
               fastq=None,
           )
       )
       wf.output("out", source=whisperalign_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for whisperAlign:

.. code-block:: bash

   # user inputs
   janis inputs whisperAlign > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fastq:
       - fastq_0.fastq.gz
       - fastq_1.fastq.gz
       index: index




5. Run whisperAlign with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       whisperAlign





Information
------------


:ID: ``whisperAlign``
:URL: `https://github.com/refresh-bio/Whisper <https://github.com/refresh-bio/Whisper>`_
:Versions: v2.0
:Container: drtomc/whisper:2.0
:Authors: Thomas Conway
:Citations: Deorowicz, S., Gudyś, A. (2019) Whisper 2: indel-sensitive short read mapping, biorXiv
:DOI: https://doi.org/10.1101/2019.12.18.881292
:Created: 2020-06-16
:Updated: 2020-06-16



Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     stdout<BAM>
======  ===========  ===============



Additional configuration (inputs)
---------------------------------

======  ===========  ========  ==========  ===========================
name    type         prefix      position  documentation
======  ===========  ========  ==========  ===========================
index   WhisperIdx                      2  base name for whisper index
fastq   FastqGzPair                     3  Paired end fastq reads
======  ===========  ========  ==========  ===========================
