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
    from janis_core.tool.tool import CommandTool, ToolArgument, ToolInput, ToolTypes, ToolOutput

    class ToolName(CommandTool):
        @staticmethod
        def tool() -> str:
            return "toolname"

        @staticmethod
        def base_command() -> Optional[Union[str, List[str]]]:
            pass

        def inputs(self) -> List[ToolInput]:
            return []

        def outputs(self) -> List[ToolOutput]:
            return []

        @staticmethod
        def container() -> str:
            pass

        @staticmethod
        def version() -> str:
            pass

        # optional

        def arguments(self) -> List[ToolArgument]:
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

- ``janis.Tool.friendly_name(self)``: A user friendly name of your tool.

- ``@staticmethod janis.Tool.tool_module()``: Unix, bioinformatics, etc.

- ``@staticmethod janis.Tool.tool_provider()``: The manafacturer of the tool, eg: Illumina, Samtools


Tool Input
**********

.. autoclass:: janis.ToolInput
   :members: __init__


Tool Output
***********

.. autoclass:: janis.ToolOutput
   :members: __init__
