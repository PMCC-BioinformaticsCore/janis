
csv
===

A comma separated file



Quickstart
-----------

.. code-block:: python

   from janis_unix.data_types.csv import Csv

   w = WorkflowBuilder("my_workflow")

   w.input("input_csv", Csv(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
