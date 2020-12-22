
FastaGzFai
==========

A local file

Secondary Files
---------------

- ``.fai``

.. note:: 

   For more information, visit `Secondary / Accessory Files <https://janis.readthedocs.io/en/latest/references/secondaryfiles.html>`__


Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.fasta import FastaGzFai

   w = WorkflowBuilder("my_workflow")

   w.input("input_fastagzfai", FastaGzFai(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
