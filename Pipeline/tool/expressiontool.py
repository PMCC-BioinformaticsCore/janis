from abc import ABC
from typing import List

from Pipeline import ToolOutput, ToolInput, String, File
from Pipeline.tool.tool import Tool, ToolTypes


class ExpressionTool(Tool, ABC):

    def id(self) -> str:
        return "ExpressionTool"

    def expression(self):
        return ""

    def cwl(self):
        pass

    def type(self):
        return ToolTypes.ExpressionTool
