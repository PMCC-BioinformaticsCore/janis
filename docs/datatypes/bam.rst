
BAM
===

A binary version of a SAM file, http://software.broadinstitute.org/software/igv/bam



Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.bam import Bam

   w = WorkflowBuilder("my_workflow")

   w.input("input_bam", Bam(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
