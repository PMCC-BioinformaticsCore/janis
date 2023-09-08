

Galaxy to CWL
=============

Currently no walkthrough tutorials for translation of Galaxy to CWL. 

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

   janis translate --from galaxy --to cwl sample_translation_files/galaxy_tool/samtools_flagstat.xml

Sample Workflow Translation

.. code-block:: bash

   janis translate --from galaxy --to cwl sample_translation_files/galaxy_workflow/rna-seq-reads-to-counts.ga


