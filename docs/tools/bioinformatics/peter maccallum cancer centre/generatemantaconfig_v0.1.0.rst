:orphan:

GenerateMantaConfig
===================

``GenerateMantaConfig`` · *1 contributor · 1 version*

Generate custom manta config file.       
        


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.generatemantaconfig import GenerateMantaConfig

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "generatemantaconfig_step",
           GenerateMantaConfig(

           )
       )
       wf.output("out", source=generatemantaconfig_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for GenerateMantaConfig:

.. code-block:: bash

   # user inputs
   janis inputs GenerateMantaConfig > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run GenerateMantaConfig with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       GenerateMantaConfig

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CodeTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          GenerateMantaConfig

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------


:ID: ``GenerateMantaConfig``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Container: python:3.8.1
:Authors: Jiaan Yu
:Citations: None
:Created: 2021-05-27 00:00:00
:Updated: 2021-05-27 00:00:00



Outputs
-----------

======  ======  ========================
name    type    documentation
======  ======  ========================
out     File    Custom Manta config file
======  ======  ========================



Additional configuration (inputs)
---------------------------------

===============  ================  ===============
name             type              documentation
===============  ================  ===============
output_filename  Optional<String>
===============  ================  ===============
    