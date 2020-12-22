
File
====

A local file



Quickstart
-----------

.. code-block:: python

   from janis_core.types.common_data_types import File

   w = WorkflowBuilder("my_workflow")

   w.input("input_file", File(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
