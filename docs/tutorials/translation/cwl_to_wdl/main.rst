

CWL to WDL 
==========

Currently no walkthrough tutorials for translation of CWL to WDL. 

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

   wget https://github.com/PMCC-BioinformaticsCore/janis/raw/translate-docs/docs/data/translation_tutorials/cwl_tool

Translate with janis

.. code-block:: bash

   janis translate --from cwl --to wdl cwl_tool/samtools_flagstat.cwl

|

Sample Workflow
---------------

Obtain source files

.. code-block:: bash

   wget https://github.com/PMCC-BioinformaticsCore/janis/raw/translate-docs/docs/data/translation_tutorials/cwl_workflow

Translate with janis

.. code-block:: bash

   janis translate --from cwl --to wdl cwl_workflow/align_sort_markdup.cwl

