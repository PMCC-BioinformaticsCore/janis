#
# Compile a java file
#

from typing import List

from Pipeline.types.common_data_types import File
from Pipeline.tool.tool import Tool, ToolInput, ToolOutput, ToolArgument


class Compile(Tool):
    file = ToolInput("file", File())
    outp = ToolOutput("outp", File(), glob="*.class")

    @staticmethod
    def tool():
        return "java-compiler"

    @staticmethod
    def base_command():
        return "javac"

    @staticmethod
    def docker():
        return "openjdk:8"

    def arguments(self) -> List[ToolArgument]:
        return [ToolArgument("$(runtime.outdir)", "-d")]
