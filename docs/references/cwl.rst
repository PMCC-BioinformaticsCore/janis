Common Workflow Language
========================

Description about the CWL project, how some of the concepts map
and limitations in the implementation.


Equivalents
-----------

- **How do I specify an Array from types?**

By wrapping the data_type by an Array, for example: ``String -> Array(String)``. Nb: the Array type must be imported from janis_core.

- **What's the equivalent for ``InitialWorkDirRequirement``?**

You can add ``localise_file=True`` to your ``ToolInput``. This is well defined for individual files in CWL and WDL. There is no equivalent for ``writable``, though suggestions are welcome. Although the ``localise_file`` attribute is allowed for Array data types, the WDL translation will become disabled as this behaviour is not well defined.

From `v1.1/CommandLineTool <https://www.commonwl.org/v1.1/CommandLineTool.html#InitialWorkDirRequirement>`_:

    *If the same File or Directory appears more than once in the InitialWorkDirRequirement listing, the implementation must choose exactly one value for path; how this value is chosen is undefined.*


- **What's the equivalent for ``EnvVarRequirement`` / how do I ensure environment variables are set in my execution environment?**

You can include the following code block within your CommandTool:

.. code-block:: python

   # Within CommandTool

   def env_vars(self):
      return {
          "myvar1": InputSelector("myInput"),
          "myvar2": "constantvalue"
      }
