
Zip
===

A zip archive, ending with .zip



Quickstart
-----------

.. code-block:: python

   from janis_unix.data_types.zipfile import ZipFile

   w = WorkflowBuilder("my_workflow")

   w.input("input_zipfile", ZipFile(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
