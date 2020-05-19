tool_template = """
from abc import ABC
from datetime import datetime
from janis_core import (
    CommandTool, ToolInput, ToolOutput, File, Boolean, 
    String, Int, Double, Float, InputSelector, Filename, 
    ToolMetadata, InputDocumentation
)

class {name}Base(CommandTool, ABC):

    def friendly_name(self) -> str:
        return "{friendly_name}"

    def tool_provider(self):
        return "{tool_provider}"

    def tool(self) -> str:
        return "{toolname}"

    def base_command(self):
        return {base_command}

    def inputs(self):
        return [
{inputs}
        ]

    def outputs(self):
        return [
{outputs}
        ]

    def metadata(self):
        return {metadata}
"""

tool_version_template = """\
from .base import {name}Base

class {name}_{escapedversion}({name}Base):
    def version(self):
        return "{version}"

    def container(self):
        return "{container}"
"""
