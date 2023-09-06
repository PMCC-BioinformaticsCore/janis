:orphan:

Parse FastQC Adapters
===========================================

``ParseFastqcAdapters`` · *2 contributors · 1 version*

Parse overrepresented region and lookup in adapter table


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.parsefastqc.v0_2_0 import ParseFastqcAdapters

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "parsefastqcadapters_step",
           ParseFastqcAdapters(
               read1_fastqc_datafile=None,
               read2_fastqc_datafile=None,
               adapters_lookup=None,
               contamination_lookup=None,
           )
       )
       wf.output("out_R1_sequences", source=parsefastqcadapters_step.out_R1_sequences)
       wf.output("out_R2_sequences", source=parsefastqcadapters_step.out_R2_sequences)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for ParseFastqcAdapters:

.. code-block:: bash

   # user inputs
   janis inputs ParseFastqcAdapters > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       adapters_lookup: adapters_lookup
       contamination_lookup: contamination_lookup
       read1_fastqc_datafile: read1_fastqc_datafile
       read2_fastqc_datafile: read2_fastqc_datafile




5. Run ParseFastqcAdapters with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       ParseFastqcAdapters

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CodeTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          ParseFastqcAdapters

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------


:ID: ``ParseFastqcAdapters``
:URL: *No URL to the documentation was provided*
:Versions: v0.2.0
:Container: python:3.8.1
:Authors: Michael Franklin, Jiaan Yu
:Citations: None
:Created: 2020-01-07 00:00:00
:Updated: 2021-10-06 00:00:00



Outputs
-----------

================  =============  ===============
name              type           documentation
================  =============  ===============
out_R1_sequences  Array<String>
out_R2_sequences  Array<String>
================  =============  ===============



Additional configuration (inputs)
---------------------------------

=====================  ======  ==========================================================================================
name                   type    documentation
=====================  ======  ==========================================================================================
read1_fastqc_datafile  File
read2_fastqc_datafile  File    fastqc_datafile of read 2
adapters_lookup        File    Specifies a file which contains the list of adapter sequences which will
                               be explicity searched against the library. The file must contain sets of named adapters in
                               the form name[tab]sequence. Lines prefixed with a hash will be ignored.
contamination_lookup   File    Specifies a file which contains the list of universal adapter
                               sequences which will be explicity searched against the library.
=====================  ======  ==========================================================================================
    