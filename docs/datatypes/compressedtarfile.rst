
CompressedTarFile
=================

A gzipped tarfile



Quickstart
-----------

.. code-block:: python

   from janis_unix.data_types.tarfile import TarFileGz

   w = WorkflowBuilder("my_workflow")

   w.input("input_tarfilegz", TarFileGz(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
