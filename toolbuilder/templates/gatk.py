gatk4_tool_template = """
from abc import ABC
from datetime import datetime
from janis_bioinformatics.tools.gatk4.gatk4toolbase import Gatk4ToolBase

from janis_core import (
    CommandTool, ToolInput, ToolOutput, File, Boolean, 
    String, Int, Double, Float, InputSelector, Filename, 
    ToolMetadata, InputDocumentation
)

class Gatk{gatk_command}Base(Gatk4ToolBase, ABC):

    @classmethod
    def gatk_command(cls):
        return "{gatk_command}"
        
    def friendly_name(self) -> str:
        return "GATK4: {gatk_command}"

    def tool(self) -> str:
        return "Gatk4{gatk_command}"

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
