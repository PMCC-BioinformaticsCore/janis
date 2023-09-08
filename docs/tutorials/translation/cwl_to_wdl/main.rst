

CWL to WDL 
==========

Currently no walkthrough tutorials for translation of CWL to WDL. 

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

   janis translate --from cwl --to wdl sample_translation_files/cwl_tool/samtools_flagstat.cwl

Sample Workflow Translation

.. code-block:: bash

   janis translate --from cwl --to wdl sample_translation_files/cwl_workflow/align_sort_markdup.cwl

