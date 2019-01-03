#
# Compile a java file
#

from typing import List

from Pipeline.types.common_data_types import File
from Pipeline.tool.commandtool import CommandTool, ToolInput, ToolOutput, ToolArgument


class Compile(CommandTool):

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

    def inputs(self):
        return [
            ToolInput("file", File())
        ]

    def outputs(self):
        return [
            ToolOutput("compiled", File(), glob="*.class")
        ]


