
CRAI
====

Index of the CRAM file https://samtools.github.io/hts-specs/CRAMv3.pdf



Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.crai import Crai

   w = WorkflowBuilder("my_workflow")

   w.input("input_crai", Crai(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
