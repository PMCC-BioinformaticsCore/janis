
jsonFile
========

A JSON file file



Quickstart
-----------

.. code-block:: python

   from janis_unix.data_types.json import JsonFile

   w = WorkflowBuilder("my_workflow")

   w.input("input_jsonfile", JsonFile(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
