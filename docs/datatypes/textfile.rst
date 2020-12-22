
TextFile
========

A textfile, ending with .txt



Quickstart
-----------

.. code-block:: python

   from janis_unix.data_types.text import TextFile

   w = WorkflowBuilder("my_workflow")

   w.input("input_textfile", TextFile(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
