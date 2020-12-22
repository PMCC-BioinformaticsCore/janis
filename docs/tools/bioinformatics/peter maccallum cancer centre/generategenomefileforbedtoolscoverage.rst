    :orphan:

    Generate genome for BedtoolsCoverage
    ============================================================================

    ``GenerateGenomeFileForBedtoolsCoverage`` · *2 contributors · 1 version*

    Generate --genome FILE for BedToolsCoverage      
        

    
    Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.generatebedtoolscoveragegenomefile import GenerateGenomeFileForBedtoolsCoverage

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "generategenomefileforbedtoolscoverage_step",
           GenerateGenomeFileForBedtoolsCoverage(
               reference=None,
           )
       )
       wf.output("out", source=generategenomefileforbedtoolscoverage_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for GenerateGenomeFileForBedtoolsCoverage:

.. code-block:: bash

   # user inputs
   janis inputs GenerateGenomeFileForBedtoolsCoverage > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run GenerateGenomeFileForBedtoolsCoverage with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       GenerateGenomeFileForBedtoolsCoverage





    Information
    ------------


    :ID: ``GenerateGenomeFileForBedtoolsCoverage``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Container: python:3.8.1
:Authors: Michael Franklin, Jiaan Yu
:Citations: None
:Created: 2020-07-21 00:00:00
:Updated: 2020-06-02 00:00:00



    Outputs
    -----------

    ======  ========  ================================
name    type      documentation
======  ========  ================================
out     TextFile  Genome file for BedToolsCoverage
======  ========  ================================



    Additional configuration (inputs)
    ---------------------------------

    ===============  ================  =====================
name             type              documentation
===============  ================  =====================
reference        FastDict
output_filename  Optional<String>  Filename to output to
===============  ================  =====================
    