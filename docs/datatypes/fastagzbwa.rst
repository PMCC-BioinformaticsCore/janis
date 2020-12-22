
FastaGzBwa
==========

A local file

Secondary Files
---------------

- ``.amb``
- ``.ann``
- ``.bwt``
- ``.pac``
- ``.sa``

.. note:: 

   For more information, visit `Secondary / Accessory Files <https://janis.readthedocs.io/en/latest/references/secondaryfiles.html>`__


Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.fasta import FastaGzBwa

   w = WorkflowBuilder("my_workflow")

   w.input("input_fastagzbwa", FastaGzBwa(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
