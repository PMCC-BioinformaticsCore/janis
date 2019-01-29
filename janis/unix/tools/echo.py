from janis import String, ToolInput, ToolOutput, Stdout
from janis.unix.tools.unixtool import UnixTool


class Echo(UnixTool):

    @staticmethod
    def tool():
        return "Echo"

    def friendly_name(self):
        return "Echo - Print to console"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [
            ToolInput("input", String())
        ]

    def outputs(self):
        return [ToolOutput("output", Stdout())]
