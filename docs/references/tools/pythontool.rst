Python Tool
============

.. note::
	Available in v0.9.0 and later

A PythonTool is a Janis mechanism for running arbitrary Python code inside a container. 

.. autoclass:: janis.PythonTool
   :members: code_block, outputs, id, version


How to build my own Tool
.........................

- Create a class that inherits from ``PythonTool``,
- The PythonTool uses the code signature for code_block to determine the inputs. This means you must use Python `type annotations`.
- The default image is ``python:3.8.1``, this is overridable by overriding the container method.

Extra notes about the code_block:

- It must be totally self contained (including ALL imports and functions)
- The only information you have available is those that you specify within your method.
- Annotation notes:
    - You can annotate with regular Python types, like: ``str`` or ``int``
    - You can annotate optionality with typings, like: ``Optional[str]`` or ``List[Optional[str]]``
    - Do NOT use ``Array(String)``, use ``List[String]`` instead.
    - Annotating with Janis types (like String, File, BamBai, etc) might cause issues with your type hints in an IDE.
- A File type will be presented to your function as a string. You will need to open the file within your function.
- You must return a dictionary with keys corresponding to the outputs of your tool.
  - File outputs should be presented as a string, and will be coerced to a File / Directory later.
- Do NOT use ``j.$TYPE`` annotations (prefixed with ``j.``, eg: ``j.File`` or ``j.String``) as annotations as this will fail at runtime.

.. code-block:: python

   from janis_core import PythonTool, TOutput, File
   from typing import Dict, Optional, List, Any
   
   class MyPythonTool(PythonTool):
       @staticmethod
       def code_block(in_string: str, in_file: File, in_integer: int) -> Dict[str, Any]:
           from shutil import copyfile
   
           copyfile(in_file, "./out.file")
               
           return {
               "myout": in_string + "-appended!",
               "myfileout": "out.file",
               "myinteger": in_integer + 1
           }
       
       def outputs(self) -> List[TOutput]:
           return [
               TOutput("myout", str),
               TOutput("myfileout", File),
               TOutput("myinteger", int)
           ]
       
       def id(self) -> str:
           return "MyPythonTool"
        
       def version(self):
           return "v0.1.0"



How to include a File inside your container
############################################

Use the type annotation ``File``, and open it yourself, for example:

.. code-block:: python

   from janis_core import PythonTool, TOutput, File
   from typing import Dict, Optional, List, Any
   
   class CatPythonTool(PythonTool):
       @staticmethod
       def code_block(my_file: File) -> Dict[str, Any]:
           
           with open(my_file, "r") as f:
               return {
                   "out": f.read()
               }

       def outputs(self) -> List[TOutput]:
           return [
               TOutput("out", str)
           ]



How to include an Array of strings inside my container
#########################################################


.. code-block:: python

   from janis_core import PythonTool, TOutput, File
   from typing import Dict, Optional, List, Any

   class JoinStrings(PythonTool):
       @staticmethod
       def code_block(my_strings: List[str], separator: str=",") -> Dict[str, Any]:

           return {"out": separator.join(my_strings)}

       def outputs(self) -> List[TOutput]:
           return [
               TOutput("out", str)
           ]