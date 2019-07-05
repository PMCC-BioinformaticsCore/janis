"""
Welcome to the source code for Janis!

Janis is a framework creating specialised, simple workflow definitions that are
transpiled to Common Workflow Language or Workflow Definition Language.

Below you'll find the classes you can access by importing Janis:

>>  import janis as j

Some noteworthy classes are Workflow, CommandTool and the janis.translations module

Some terminology:
    - Edge:
        - An edge may exist between 2 nodes, it represents a dependency
        - Every edge will have a source_map which is indexed by the tag of the input on the finish node
        - This source_map value is singular or a list, as you can connect multiple sources to one input



_Proudly made on Planet Earth._

"""

# PEP396:  https://www.python.org/dev/peps/pep-0396/
from janis.__meta__ import __version__

from janis.hints import CaptureType, Engine, HINTS, Hint, HintEnum, HintArray
from janis.tool.commandtool import CommandTool
from janis.tool.tool import Tool, ToolArgument, ToolInput, ToolOutput
from janis.translations import SupportedTranslations
from janis.types import (
    InputSelector,
    WildcardSelector,
    MemorySelector,
    CpuSelector,
    StringFormatter,
)
from janis.types.common_data_types import (
    Boolean,
    String,
    Int,
    Float,
    Double,
    File,
    Directory,
    Array,
    Filename,
    Stdout,
)
from janis.types.data_types import DataType
from janis.unix import *
from janis.utils.logger import Logger, LogLevel
from janis.utils.metadata import Metadata, WorkflowMetadata, ToolMetadata
from janis.workflow.input import Input
from janis.workflow.output import Output
from janis.workflow.step import Step
from janis.workflow.workflow import Workflow
