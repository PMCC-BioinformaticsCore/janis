Code Tool
==========

.. note::
	BETA: Available in v0.9.0 and later

A code tool is an abstract concept in Janis, that aims to execute arbitrary code inside a container.

Currently there is only one type of CodeTool:

- :class:`janis.PythonTool`

The aim is to make it simpler to perform basic steps within your container, for example say I want to perform some basic processing of a file that isn't trivial in Janis (and hence CWL / WDL) but is in a programming language, Janis gives you the functionality to make this possible.


Creating a new code tool type
.............................

This process is designed to be fairly simple, but there are a few important notes:

- Your tool must log ONLY a JSON string to stdout (this will get parsed later)
- The ``prepared_script`` block must return a string that will get written to a file within the container to be executed that can accept inputs via command line. Within the PythonTool, you'll see that a parser (by ``argparse``) is generated in this method.


.. code-block:: python

   import janis_core as j

   class LanguageTool(CodeTool):
       # You might leave these fields to be overriden by the user
       def inputs(self) -> List[TInput]:
           pass
   
       def outputs(self) -> List[TOutput]:
           pass
   
       def base_command(self):
           pass
   
       def script_name(self):
           pass
   
       def container(self):
           pass
   
       def prepared_script(self):
           pass
   
       def id(self) -> str:
           pass
   
       def version(self):
           pass

