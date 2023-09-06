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

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CodeTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          prepareALLSortsInput

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











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
    