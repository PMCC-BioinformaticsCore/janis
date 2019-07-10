#
# Compile a java file
#

from typing import List

from janis import ToolInput, ToolOutput, ToolArgument, WildcardSelector, File
from janis.unix.tools.unixtool import UnixTool


class Compile(UnixTool):
    @staticmethod
    def tool():
        return "javaCompiler"

    def friendly_name(self):
        return "Java compiler"

    @staticmethod
    def base_command():
        return "javac"

    @staticmethod
    def docker():
        return "openjdk:8"

    def arguments(self) -> List[ToolArgument]:
        return [ToolArgument(".", "-d")]  # CurrentWorkingDirectory()

    def inputs(self):
        return [ToolInput("file", File(), position=1)]

    def outputs(self):
        return [ToolOutput("out", File(), glob=WildcardSelector("*.class"))]
