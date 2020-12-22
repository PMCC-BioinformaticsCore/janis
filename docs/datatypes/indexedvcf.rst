
IndexedVCF
==========

Variant Call Format:

    The Variant Call Format (VCF) specifies the format of a text file 
    used in bioinformatics for storing gene sequence variations. 

    Documentation: https://samtools.github.io/hts-specs/VCFv4.3.pdf

Secondary Files
---------------

- ``.idx``

.. note:: 

   For more information, visit `Secondary / Accessory Files <https://janis.readthedocs.io/en/latest/references/secondaryfiles.html>`__


Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.vcf import VcfIdx

   w = WorkflowBuilder("my_workflow")

   w.input("input_vcfidx", VcfIdx(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
