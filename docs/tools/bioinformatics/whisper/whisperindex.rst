:orphan:

Whisper-Index
============================

*1 contributor · 1 version*

Builds a whisper index


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.whisper.index.versions import WhisperIndex_2_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "whisperindex_step",
           WhisperIndex_2_0(
               index_name=None,
               fasta=None,
           )
       )
       wf.output("out", source=whisperindex_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for WhisperIndex:

.. code-block:: bash

   # user inputs
   janis inputs WhisperIndex > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fasta:
       - fasta_0.fasta
       - fasta_1.fasta
       index_name: <value>




5. Run WhisperIndex with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       WhisperIndex





Information
------------


:ID: ``WhisperIndex``
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

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     WhisperIdx
======  ==========  ===============



Additional configuration (inputs)
---------------------------------

==========  ============  ========  ==========  ====================
name        type          prefix      position  documentation
==========  ============  ========  ==========  ====================
index_name  String                           2  name of the index
fasta       Array<Fasta>                     3  FASTA files to index
==========  ============  ========  ==========  ====================
