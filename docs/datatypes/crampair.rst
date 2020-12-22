
CramPair
========

A Cram and Crai as the secondary

Secondary Files
---------------

- ``.crai``

.. note:: 

   For more information, visit `Secondary / Accessory Files <https://janis.readthedocs.io/en/latest/references/secondaryfiles.html>`__


Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.cram import CramCrai

   w = WorkflowBuilder("my_workflow")

   w.input("input_cramcrai", CramCrai(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
