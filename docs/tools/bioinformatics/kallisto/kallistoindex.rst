:orphan:

Kallisto-Index
==============================

*1 contributor · 1 version*

Builds a kallisto index


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.kallisto.index.versions import KallistoIndex_0_46_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "kallistoindex_step",
           KallistoIndex_0_46_2(
               reference=None,
           )
       )
       wf.output("out", source=kallistoindex_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for kallistoIndex:

.. code-block:: bash

   # user inputs
   janis inputs kallistoIndex > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run kallistoIndex with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       kallistoIndex





Information
------------


:ID: ``kallistoIndex``
:URL: `https://pachterlab.github.io/kallisto/manual.html <https://pachterlab.github.io/kallisto/manual.html>`_
:Versions: v0.46.2
:Container: quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1
:Authors: Thomas Conway
:Citations: NL Bray, H Pimentel, P Melsted and L Pachter, Near optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, p 525--527 (2016).
:DOI: https://doi.org/10.1038/nbt.3519
:Created: 2020-05-25
:Updated: 2020-05-25



Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     KallistoIdx
======  ===========  ===============



Additional configuration (inputs)
---------------------------------

=========  ==================  ========  ==========  =================================================
name       type                prefix      position  documentation
=========  ==================  ========  ==========  =================================================
reference  Fasta                                  3  Filename for a reference transcriptome
kmer_size  Optional<Integer>   -k                 1  k-mer (odd) length (default: 31, max value: 31)
index      Optional<Filename>  -i                 2  Filename for the kallisto index to be constructed
=========  ==================  ========  ==========  =================================================
