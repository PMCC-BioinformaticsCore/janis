
FastDict
========

A local file

Secondary Files
---------------

- ``^.dict``

.. note:: 

   For more information, visit `Secondary / Accessory Files <https://janis.readthedocs.io/en/latest/references/secondaryfiles.html>`__


Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.fasta import FastaDict

   w = WorkflowBuilder("my_workflow")

   w.input("input_fastadict", FastaDict(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
