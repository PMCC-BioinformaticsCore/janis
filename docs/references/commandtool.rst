Command Tool
=============

*A class that wraps a CommandLineTool with named argument*

Declaration
***********

.. autoclass:: janis.CommandTool

Overview
*********

A :class:`janis.CommandTool` is the programmatic interface between a
set of inputs, how these inputs bind on the Command Line to call the
tool, and how to collect the outputs. It would be rare that you should
directly instantiate ``CommandTool``, instead you should subclass and
and override the methods you need as declared below.


Template
.........

.. code-block:: python

   from typing import List, Optional, Union
   import janis_core as j

   class ToolName(j.CommandTool):
       @staticmethod
       def tool() -> str:
           return "toolname"

       @staticmethod
       def base_command() -> Optional[Union[str, List[str]]]:
           pass

       def inputs(self) -> List[j.ToolInput]:
           return []

       def outputs(self) -> List[j.ToolOutput]:
           return []

       @staticmethod
       def container() -> str:
           return ""

       @staticmethod
       def version() -> str:
           pass

        # optional

        def arguments(self) -> List[j.ToolArgument]:
            return []

        def friendly_name(self) -> str:
            pass

        @staticmethod
        def tool_module() -> str:
            pass

        @staticmethod
        def tool_provider() -> str:
            pass

Structure
..........

A new tool definition must subclass the `janis.CommandTool` class and implement the required abstract methods:

- ``@staticmethod janis.CommandTool.tool() -> str``: Unique identifier of the tool

- ``@staticmethod janis.CommandTool.base_command() -> str``: The command of the tool to execute, usually the tool name or path and not related to any inputs.

- ``janis.CommandTool.inputs(self) -> List[``:class:`janis.ToolInput` ``]``: A list of named tool inputs that will be used to create the command line.

- ``janis.CommandTool.arguments(self) -> List[``:class:`janis.ToolArgument` ``]``: A list of arguments that will be used to create the command line.

- ``janis.CommandTool.outputs(self) -> List[``:class:`janis.ToolOutput` ``]``: A list of named outputs of the tool; a ``ToolOutput`` declares how the output is captured.

- ``@staticmethod janis.CommandTool.container()``: A link to an OCI compliant container accessible by the engine. Previously, ``docker()``.

- ``@staticmethod janis.CommandTool.version()``: Version of the tool.

To better categorise your tool, you can additionally implement the following methods:

- ``janis.Tool.friendly_name(self)``: A user friendly name of your tool (must be implemented for generated docs)

- ``@staticmethod janis.Tool.tool_module()``: Unix, bioinformatics, etc.

- ``@staticmethod janis.Tool.tool_provider()``: The manafacturer of the tool, eg: Illumina, Samtools


Tool Input
**********

.. autoclass:: janis.ToolInput
   :members: __init__

.. note::

   A ToolInput (and ``ToolArgument``) must have either a ``position`` AND / OR ``prefix`` to be bound onto the command line.

- The ``prefix`` is a string that precedes the inputted value. By default the prefix is separated by a space, however this can be removed with ``separate_value_from_prefix=False``.
- The ``position`` represents the order of how arguments are bound onto the command line. Lower numbers get a higher priority, not providing a number will default to 0.
- ``prefix_applies_to_all_elements`` applies the prefix to each element in an array (only applicable for array inputs).


Tool Output
***********

.. autoclass:: janis.ToolOutput
   :members: __init__
