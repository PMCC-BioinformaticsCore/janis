
Gzip
====

A gzipped file



Quickstart
-----------

.. code-block:: python

   from janis_unix.data_types.gunzipped import Gunzipped

   w = WorkflowBuilder("my_workflow")

   w.input("input_gunzipped", Gunzipped(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
