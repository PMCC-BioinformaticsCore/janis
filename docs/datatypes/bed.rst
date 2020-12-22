
bed
===

A local file



Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.bed import Bed

   w = WorkflowBuilder("my_workflow")

   w.input("input_bed", Bed(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
