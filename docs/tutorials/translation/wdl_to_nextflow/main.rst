

WDL to Nextflow
===============

Currently no walkthrough tutorials for translation of WDL to Nextflow. 

Run the following to see janis in action:

Install Janis

.. code-block:: bash

   python -m venv venv 
   source venv/bin/activate
   pip install janis-pipelines

Obtain Sample Files

.. code-block:: bash

   git clone https://github.com/GraceAHall/sample_translation_files

Sample Tool Translation

.. code-block:: bash

   janis translate --from wdl --to nextflow sample_translation_files/wdl_tool/cutadapt.wdl

Sample Workflow Translation

.. code-block:: bash

   janis translate --from wdl --to nextflow sample_translation_files/wdl_workflow/linear.wdl

