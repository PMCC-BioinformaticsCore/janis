from janis import String, ToolInput, ToolOutput, Stdout
from janis.unix.tools.unixtool import UnixTool


class Echo(UnixTool):
    @staticmethod
    def tool():
        return "echo"

    def friendly_name(self):
        return "Echo - Print to console"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [ToolInput("inp", String(), position=0)]

    def outputs(self):
        return [ToolOutput("out", Stdout())]
