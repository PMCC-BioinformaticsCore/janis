
FastaGz
=======

A local file



Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.fasta import FastaGz

   w = WorkflowBuilder("my_workflow")

   w.input("input_fastagz", FastaGz(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
