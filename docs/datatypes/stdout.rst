
Stdout
======

A local file



Quickstart
-----------

.. code-block:: python

   from janis_core.types.common_data_types import Stdout

   w = WorkflowBuilder("my_workflow")

   w.input("input_stdout", Stdout(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
