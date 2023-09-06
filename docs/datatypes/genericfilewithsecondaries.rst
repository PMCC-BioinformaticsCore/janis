
GenericFileWithSecondaries
==========================

A local file

Secondary Files
---------------

- ``None``

.. note:: 

   For more information, visit `Secondary / Accessory Files <https://janis.readthedocs.io/en/latest/references/secondaryfiles.html>`__


Quickstart
-----------

.. code-block:: python

   from janis_core.types.common_data_types import GenericFileWithSecondaries

   w = WorkflowBuilder("my_workflow")

   w.input("input_genericfilewithsecondaries", GenericFileWithSecondaries(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2023-09-06*.
