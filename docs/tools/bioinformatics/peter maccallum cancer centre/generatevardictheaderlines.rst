    :orphan:

    GenerateVardictHeaderLines
    ==========================

    ``GenerateVardictHeaderLines`` · *2 contributors · 1 version*

    Generate VarDict Headerlines.       
        

    
    Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.generatevardictheaderlines import GenerateVardictHeaderLines

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "generatevardictheaderlines_step",
           GenerateVardictHeaderLines(
               reference=None,
           )
       )
       wf.output("out", source=generatevardictheaderlines_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for GenerateVardictHeaderLines:

.. code-block:: bash

   # user inputs
   janis inputs GenerateVardictHeaderLines > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run GenerateVardictHeaderLines with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       GenerateVardictHeaderLines





    Information
    ------------


    :ID: ``GenerateVardictHeaderLines``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Container: python:3.8.1
:Authors: Michael Franklin, Jiaan Yu
:Citations: None
:Created: 2020-06-02 00:00:00
:Updated: 2020-06-02 00:00:00



    Outputs
    -----------

    ======  ======  ===============================================================
name    type    documentation
======  ======  ===============================================================
out     File    Header file for VarDict, generated based on the reference index
======  ======  ===============================================================



    Additional configuration (inputs)
    ---------------------------------

    ===============  ================  =====================
name             type              documentation
===============  ================  =====================
reference        FastDict
output_filename  Optional<String>  Filename to output to
===============  ================  =====================
    