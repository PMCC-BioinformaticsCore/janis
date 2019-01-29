#
# Compile a java file
#

from typing import List

from Pipeline.types.common_data_types import File, CurrentWorkingDirectory
from Pipeline.tool.commandtool import CommandTool, ToolInput, ToolOutput, ToolArgument


class Compile(CommandTool):

    @staticmethod
    def tool():
        return "javaCompiler"

    @staticmethod
    def base_command():
        return "javac"

    def friendly_name(self):
        return "Java compiler"

    @staticmethod
    def docker():
        return "openjdk:8"

    def arguments(self) -> List[ToolArgument]:
        return [ToolArgument(".", "-d")]    # CurrentWorkingDirectory()

    def inputs(self):
        return [
            ToolInput("file", File(), position=1)
        ]

    def outputs(self):
        return [
            ToolOutput("compiled", File(), glob="*.class")
        ]


