
tsv
===

A tab separated file



Quickstart
-----------

.. code-block:: python

   from janis_unix.data_types.tsv import Tsv

   w = WorkflowBuilder("my_workflow")

   w.input("input_tsv", Tsv(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
