    :orphan:

    Create genomic call regions
    ===============================================

    ``CreateCallRegions`` · *1 contributor · 1 version*

    No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

    
    Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.dawson.createcallregions.base import CreateCallRegions

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "createcallregions_step",
           CreateCallRegions(
               reference=None,
               regionSize=None,
           )
       )
       wf.output("regions", source=createcallregions_step.regions)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for CreateCallRegions:

.. code-block:: bash

   # user inputs
   janis inputs CreateCallRegions > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta
       regionSize: 0




5. Run CreateCallRegions with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       CreateCallRegions





    Information
    ------------


    :ID: ``CreateCallRegions``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Container: python:3.8.1
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2020-06-17
:Updated: 2020-07-16



    Outputs
    -----------

    =======  =============  ===============
name     type           documentation
=======  =============  ===============
regions  Array<String>
=======  =============  ===============



    Additional configuration (inputs)
    ---------------------------------

    ==========  =================  ===============
name        type               documentation
==========  =================  ===============
reference   FastaFai
regionSize  Integer
equalize    Optional<Boolean>
==========  =================  ===============
    