Command Tool
=============

*A class that wraps a CommandLineTool with named argument*

Declaration
***********

.. autoclass:: janis.CommandTool

.. autoclass:: janis.CommandToolBuilder
   :members: __init__


Overview
*********

A :class:`janis.CommandTool` is the programmatic interface between a
set of inputs, how these inputs bind on the Command Line to call the
tool, and how to collect the outputs. It would be rare that you should
directly instantiate ``CommandTool``, instead you should subclass and
and override the methods you need as declared below.

Like a workflow, there are two methods to declare a command tool:

- Use the ``CommandToolBuilder`` class
- Inherit from the ``CommandTool`` class,


Template
.........

CommandToolBuilder:

.. code-block:: python

   import janis_core as j

   ToolName = j.CommandToolBuilder(
       tool: str="toolname",
       base_command=["base", "command"],
       inputs: List[j.ToolInput]=[],
       outputs: List[j.ToolOutput]=[],
       container="container/name:version",
       version="version",
       friendly_name=None,
       arguments=None,
       env_vars=None,
       tool_module=None,
       tool_provider=None,
       metadata: ToolMetadata=j.ToolMetadata(),
       cpu: Union[int, Callable[[Dict[str, Any]], int]]=None,
       memory: Union[int, Callable[[Dict[str, Any]], int]]=None,
   )


This is equivalent to the inherited template:

.. code-block:: python

   from typing import List, Optional, Union
   import janis_core as j

   class ToolName(j.CommandTool):
       def tool(self) -> str:
           return "toolname"

       def base_command(self) -> Optional[Union[str, List[str]]]:
           pass

       def inputs(self) -> List[j.ToolInput]:
           return []

       def outputs(self) -> List[j.ToolOutput]:
           return []

       def container(self) -> str:
           return ""

       def version(self) -> str:
           pass

        # optional

        def arguments(self) -> List[j.ToolArgument]:
            return []

        def env_vars(self) -> Optional[Dict[str, Union[str, Selector]]]:
            return {}

        def friendly_name(self) -> str:
            pass

        def tool_module(self) -> str:
            pass

        def tool_provider(self) -> str:
            pass

        def cpu(self, hints: Dict) -> int:
            pass

        def memory(self, hints: Dict) -> int:
            pass

        def bind_metadata(self) -> j.ToolMetadata:
            pass

Structure
..........

A new tool definition must subclass the `janis.CommandTool` class and implement the required abstract methods:

- ``janis.CommandTool.tool(self) -> str``: Unique identifier of the tool

- ``janis.CommandTool.base_command(self) -> str``: The command of the tool to execute, usually the tool name or path and not related to any inputs.

- ``janis.CommandTool.inputs(self) -> List[``:class:`janis.ToolInput` ``]``: A list of named tool inputs that will be used to create the command line.

- ``janis.CommandTool.arguments(self) -> List[``:class:`janis.ToolArgument` ``]``: A list of arguments that will be used to create the command line.

- ``janis.CommandTool.outputs(self) -> List[``:class:`janis.ToolOutput` ``]``: A list of named outputs of the tool; a ``ToolOutput`` declares how the output is captured.

- ``janis.CommandTool.container(self)``: A link to an OCI compliant container accessible by the engine. Previously, ``docker()``.

- ``janis.CommandTool.version(self)``: Version of the tool.

- ``janis.CommandTool.env_vars(self) -> Optional[Dict[str, Union[str, Selector]]]``: A dictionary of environment variables that should be defined within the container.

To better categorise your tool, you can additionally implement the following methods:

- ``janis.Tool.friendly_name(self)``: A user friendly name of your tool (must be implemented for generated docs)

- ``janis.Tool.tool_module(self)``: Unix, bioinformatics, etc.

- ``janis.Tool.tool_provider(self)``: The manafacturer of the tool, eg: Illumina, Samtools


Tool Input
**********

.. autoclass:: janis.ToolInput
   :members: __init__

.. note::

   A ToolInput (and ``ToolArgument``) must have either a ``position`` AND / OR ``prefix`` to be bound onto the command line.

- The ``prefix`` is a string that precedes the inputted value. By default the prefix is separated by a space, however this can be removed with ``separate_value_from_prefix=False``.
- The ``position`` represents the order of how arguments are bound onto the command line. Lower numbers get a higher priority, not providing a number will default to 0.
- ``prefix_applies_to_all_elements`` applies the prefix to each element in an array (only applicable for array inputs).
- The ``localise_file`` attribute places the file input within the execution directory.
- ``presents_as`` is a mechanism for overriding the name to localise to. The ``localise_file`` parameter MUST be set to `True` for ``presents_as``
- ``secondaries_present_as`` is a mechanism for overriding the format of secondary files. ``localise_file`` does NOT need to be set for this functionality to work. In CWL, this relies on https://github.com/common-workflow-language/cwltool/pull/1233

Tool Argument
**************

.. autoclass:: janis.ToolArgument
   :members: __init__


Tool Output
***********

.. autoclass:: janis.ToolOutput
   :members: __init__
