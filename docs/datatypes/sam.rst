
SAM
===

Tab-delimited text file that contains sequence alignment data



Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.sam import Sam

   w = WorkflowBuilder("my_workflow")

   w.input("input_sam", Sam(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
