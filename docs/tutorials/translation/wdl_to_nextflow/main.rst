

WDL to Nextflow
===============

Currently no walkthrough tutorials for translation of WDL to Nextflow. 

Run the following to see janis in action:

Installing Janis

.. code-block:: bash

   python -m venv venv 
   source venv/bin/activate
   pip install janis-pipelines

|

Sample Tool
-----------

Obtain source workflow files

.. code-block:: bash

   wget https://github.com/PMCC-BioinformaticsCore/janis/raw/translate-docs/docs/data/translation_tutorials/wdl_tool/cutadapt.wdl

Translate with janis

.. code-block:: bash

   janis translate --from wdl --to nextflow cutadapt.wdl

|

Sample Workflow
---------------

Obtain source workflow files

.. code-block:: bash

   wget https://github.com/PMCC-BioinformaticsCore/janis/raw/translate-docs/docs/data/translation_tutorials/wdl_workflow/linear.wdl

Translate with janis

.. code-block:: bash

   janis translate --from wdl --to nextflow linear.wdl

