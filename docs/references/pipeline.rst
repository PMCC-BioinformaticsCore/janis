Pipeline
********

The pipeline module contains all the classes you should need to build a workflow. Technically it also contains all the tool definitions you would need but they're nested under different git submodules.

These classes are based on a superset of the features available within CWL and WDL. You should have a fairly good understanding about workflow concepts such as: Inputs, Steps (+ tools) and Outputs.

If you want to build tool definitions, you will need a moderate understanding of the command line, how parameters are organised and knowledge about the tool.


Building workflows
==================

.. autoclass:: Pipeline.Workflow
.. autoclass:: Pipeline.Input
.. autoclass:: Pipeline.Output
.. autoclass:: Pipeline.Step


Data types
==========

.. autoclass:: Pipeline.String
.. autoclass:: Pipeline.File
.. autoclass:: Pipeline.Int
.. autoclass:: Pipeline.Float
.. autoclass:: Pipeline.File
