
Stderr
======

A local file



Quickstart
-----------

.. code-block:: python

   from janis_core.types.common_data_types import Stderr

   w = WorkflowBuilder("my_workflow")

   w.input("input_stderr", Stderr(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-11-10*.
