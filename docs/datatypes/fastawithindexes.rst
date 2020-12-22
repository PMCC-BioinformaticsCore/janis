
FastaWithIndexes
================

A local file

Secondary Files
---------------

- ``.fai``
- ``.amb``
- ``.ann``
- ``.bwt``
- ``.pac``
- ``.sa``
- ``^.dict``

.. note:: 

   For more information, visit `Secondary / Accessory Files <https://janis.readthedocs.io/en/latest/references/secondaryfiles.html>`__


Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.fasta import FastaWithIndexes

   w = WorkflowBuilder("my_workflow")

   w.input("input_fastawithindexes", FastaWithIndexes(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
