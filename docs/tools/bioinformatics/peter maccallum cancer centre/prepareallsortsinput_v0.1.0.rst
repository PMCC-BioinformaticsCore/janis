    :orphan:

    Prepare ALLSorts Input
    =============================================

    ``prepareALLSortsInput`` · *1 contributor · 1 version*

    No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

    
    Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.oshlack.prepareallsortsinput import PrepareALLSortsInput_0_1_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "prepareallsortsinput_step",
           PrepareALLSortsInput_0_1_0(
               inps=None,
           )
       )
       wf.output("out", source=prepareallsortsinput_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for prepareALLSortsInput:

.. code-block:: bash

   # user inputs
   janis inputs prepareALLSortsInput > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inps:
       - inps_0
       - inps_1




5. Run prepareALLSortsInput with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       prepareALLSortsInput





    Information
    ------------


    :ID: ``prepareALLSortsInput``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Container: python:3.8.1
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-09-21
:Updated: 2020-09-21



    Outputs
    -----------

    ======  ======  ===============
name    type    documentation
======  ======  ===============
out     csv
======  ======  ===============



    Additional configuration (inputs)
    ---------------------------------

    ===============  =======================  ===============
name             type                     documentation
===============  =======================  ===============
inps             Array<File>
labels           Optional<Array<String>>
output_filename  Optional<String>
fusion_caller    Optional<String>
===============  =======================  ===============
    