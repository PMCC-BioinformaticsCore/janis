
CRAM
====

A binary version of a SAM file, https://samtools.github.io/hts-specs/CRAMv3.pdf



Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.cram import Cram

   w = WorkflowBuilder("my_workflow")

   w.input("input_cram", Cram(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
