    :orphan:

    Generating genomic intervals by chromosome
    ==========================================================================

    ``GenerateIntervalsByChromosome`` · *1 contributor · 1 version*

    No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

    
    Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.generatintervalsbychromosome.generateintervalsbychromosome import GenerateIntervalsByChromosome

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "generateintervalsbychromosome_step",
           GenerateIntervalsByChromosome(
               reference=None,
           )
       )
       wf.output("out_regions", source=generateintervalsbychromosome_step.out_regions)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for GenerateIntervalsByChromosome:

.. code-block:: bash

   # user inputs
   janis inputs GenerateIntervalsByChromosome > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run GenerateIntervalsByChromosome with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       GenerateIntervalsByChromosome





    Information
    ------------


    :ID: ``GenerateIntervalsByChromosome``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Container: python:3.8.1
:Authors: Michael Franklin
:Citations: None
:Created: 2020-10-19 00:00:00
:Updated: 2020-10-19 00:00:00



    Outputs
    -----------

    ===========  ==========  ===============
name         type        documentation
===========  ==========  ===============
out_regions  Array<bed>
===========  ==========  ===============



    Additional configuration (inputs)
    ---------------------------------

    ===============  =======================  =========================================================================
name             type                     documentation
===============  =======================  =========================================================================
reference        FastDict                 FASTA reference with ^.dict reference
prefix           Optional<String>         contig prefix, default 'chr'
allowed_contigs  Optional<Array<String>>  Limits allowed_contigs to a list, this defaults of Human CHRs, 1-23,X,Y,Z
max_size         Optional<Integer>        Max size of interval, maybe 5000 for VarDict.
overlap          Optional<Integer>        Consider indels spanning regions, so choose
single_file      Optional<Boolean>        Produce a SINGLE .bed file with all the listed regions
===============  =======================  =========================================================================
    