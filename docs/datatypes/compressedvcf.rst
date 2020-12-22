
CompressedVCF
=============

.vcf.gz



Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.vcf import CompressedVcf

   w = WorkflowBuilder("my_workflow")

   w.input("input_compressedvcf", CompressedVcf(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
