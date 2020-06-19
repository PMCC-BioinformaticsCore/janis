:orphan:

Kallisto-Quant
==============================

*1 contributor · 1 version*

Builds a kallisto index


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.kallisto.quant.versions import KallistoQuant_0_46_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "kallistoquant_step",
           KallistoQuant_0_46_2(
               index=None,
               fastq=None,
           )
       )
       wf.output("out", source=kallistoquant_step.out)
       wf.output("stats", source=kallistoquant_step.stats)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for kallistoQuant:

.. code-block:: bash

   # user inputs
   janis inputs kallistoQuant > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fastq:
       - fastq_0.fastq
       - fastq_1.fastq
       index: index.kidx




5. Run kallistoQuant with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       kallistoQuant





Information
------------


:ID: ``kallistoQuant``
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

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
stats   File
======  ======  ===============



Additional configuration (inputs)
---------------------------------

===============  ==================  =================  ==========  ========================================================================================
name             type                prefix               position  documentation
===============  ==================  =================  ==========  ========================================================================================
index            KallistoIdx         -i                          2  Filename for the kallisto index to be constructed
fastq            Array<Fastq>                                    4  FASTQ files to process
outdir           Optional<Filename>  -o                          3  directory to put outputs in
bias             Optional<Boolean>   --bias                         Perform sequence based bias correction
fusion           Optional<Boolean>   --fusion                       Search for fusions for Pizzly
single           Optional<Boolean>   --single                       Quantify single-end reads
overhang         Optional<Boolean>   --single-overhang              Include reads where unobserved rest of fragment is predicted to lie outside a transcript
fr_stranded      Optional<Boolean>   --fr-stranded                  Strand specific reads, first read forward
rf_stranded      Optional<Boolean>   --rf-stranded                  Strand specific reads, first read reverse
fragment_length  Optional<Double>    -l                             Estimated average fragment length
fragment_sd      Optional<Double>    -s                             Estimated standard deviation of fragment length
===============  ==================  =================  ==========  ========================================================================================
