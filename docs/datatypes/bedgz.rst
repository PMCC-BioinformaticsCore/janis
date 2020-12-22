
BedGz
=====

A gzipped file



Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.bed import BedGz

   w = WorkflowBuilder("my_workflow")

   w.input("input_bedgz", BedGz(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
