
TarFile
=======

A tarfile, ending with .tar



Quickstart
-----------

.. code-block:: python

   from janis_unix.data_types.tarfile import TarFile

   w = WorkflowBuilder("my_workflow")

   w.input("input_tarfile", TarFile(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
