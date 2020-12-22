
IndexedBam
==========

A Bam and bai as the secondary

Secondary Files
---------------

- ``.bai``

.. note:: 

   For more information, visit `Secondary / Accessory Files <https://janis.readthedocs.io/en/latest/references/secondaryfiles.html>`__


Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.bam import BamBai

   w = WorkflowBuilder("my_workflow")

   w.input("input_bambai", BamBai(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
