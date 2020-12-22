
BedTABIX
========

A gzipped file

Secondary Files
---------------

- ``.tbi``

.. note:: 

   For more information, visit `Secondary / Accessory Files <https://janis.readthedocs.io/en/latest/references/secondaryfiles.html>`__


Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.bed import BedTabix

   w = WorkflowBuilder("my_workflow")

   w.input("input_bedtabix", BedTabix(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
