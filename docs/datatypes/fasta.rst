
Fasta
=====

A local file



Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.fasta import Fasta

   w = WorkflowBuilder("my_workflow")

   w.input("input_fasta", Fasta(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
