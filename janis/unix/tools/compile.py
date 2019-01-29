#
# Compile a java file
#

from typing import List

from janis.tool.commandtool import ToolInput, ToolOutput, ToolArgument
from janis.types.common_data_types import File
from janis.unix.tools.unixtool import UnixTool


class Compile(UnixTool):

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


