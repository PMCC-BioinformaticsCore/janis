    :orphan:

    Parse FastQC Adaptors
    ===========================================

    ``ParseFastqcAdaptors`` · *1 contributor · 1 version*

    Parse overrepresented region and lookup in Cutadapt table

    
    Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.parsefastqc.v0_1_0 import ParseFastqcAdaptors

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "parsefastqcadaptors_step",
           ParseFastqcAdaptors(
               fastqc_datafiles=None,
           )
       )
       wf.output("adaptor_sequences", source=parsefastqcadaptors_step.adaptor_sequences)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for ParseFastqcAdaptors:

.. code-block:: bash

   # user inputs
   janis inputs ParseFastqcAdaptors > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fastqc_datafiles:
       - fastqc_datafiles_0
       - fastqc_datafiles_1




5. Run ParseFastqcAdaptors with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       ParseFastqcAdaptors





    Information
    ------------


    :ID: ``ParseFastqcAdaptors``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Container: python:3.8.1
:Authors: Michael Franklin
:Citations: None
:Created: 2020-01-07 00:00:00
:Updated: 2020-02-14 00:00:00



    Outputs
    -----------

    =================  =============  ===============
name               type           documentation
=================  =============  ===============
adaptor_sequences  Array<String>
=================  =============  ===============



    Additional configuration (inputs)
    ---------------------------------

    ========================  ==============  ==========================================================================================
name                      type            documentation
========================  ==============  ==========================================================================================
fastqc_datafiles          Array<File>
cutadapt_adaptors_lookup  Optional<File>  Specifies a file which contains the list of adapter sequences which will
                                          be explicity searched against the library. The file must contain sets of named adapters in
                                          the form name[tab]sequence. Lines prefixed with a hash will be ignored.
========================  ==============  ==========================================================================================
    