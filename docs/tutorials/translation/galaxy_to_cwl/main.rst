

Galaxy to CWL
=============

Currently no walkthrough tutorials for translation of Galaxy to CWL. 

Run the following to see janis in action:

Installation

.. code-block:: bash

   python -m venv venv 
   source venv/bin/activate
   pip install janis-pipelines

|

Sample Tool
-----------

Obtain source files

.. code-block:: bash

   wget https://github.com/PMCC-BioinformaticsCore/janis/raw/translate-docs/docs/data/translation_tutorials/galaxy_tool

Translate with janis

.. code-block:: bash

   janis translate --from galaxy --to cwl galaxy_tool/samtools_flagstat.xml

|

Sample Workflow
---------------

Obtain source files

.. code-block:: bash

   wget https://github.com/PMCC-BioinformaticsCore/janis/raw/translate-docs/docs/data/translation_tutorials/galaxy_workflow

Translate with janis

.. code-block:: bash

   janis translate --from galaxy --to cwl galaxy_workflow/rna-seq-reads-to-counts.ga

