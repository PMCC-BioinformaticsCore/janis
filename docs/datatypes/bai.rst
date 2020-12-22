
BAI
===

Index of the BAM file (https://www.biostars.org/p/15847/), http://software.broadinstitute.org/software/igv/bam



Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.bai import Bai

   w = WorkflowBuilder("my_workflow")

   w.input("input_bai", Bai(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
