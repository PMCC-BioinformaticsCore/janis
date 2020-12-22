
CompressedIndexedVCF
====================

.vcf.gz with .vcf.gz.tbi file

Secondary Files
---------------

- ``.tbi``

.. note:: 

   For more information, visit `Secondary / Accessory Files <https://janis.readthedocs.io/en/latest/references/secondaryfiles.html>`__


Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.vcf import VcfTabix

   w = WorkflowBuilder("my_workflow")

   w.input("input_vcftabix", VcfTabix(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
