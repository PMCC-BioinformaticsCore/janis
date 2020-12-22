
Filename
========

This class is a placeholder for generated filenames, by default it is optional and CAN be overridden, 
however the program has been structured in a way such that these names will be generated based on the step label. 
These should only be used when the tool _requires_ a filename to output and you aren't 
concerned what the filename should be. The Filename DataType should NOT be used as an output.



Quickstart
-----------

.. code-block:: python

   from janis_core.types.common_data_types import Filename

   w = WorkflowBuilder("my_workflow")

   w.input("input_filename", Filename(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
