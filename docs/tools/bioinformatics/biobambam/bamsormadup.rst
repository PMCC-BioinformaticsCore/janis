:orphan:

BamSorMaDup
=========================

*1 contributor · 1 version*

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
