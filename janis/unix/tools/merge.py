from janis import File, Array, ToolInput, ToolOutput, Stdout
from janis.unix.tools.unixtool import UnixTool


class Merge(UnixTool):

    def friendly_name(self) -> str:
        return "Merge Files"

    @staticmethod
    def tool():
        return "merge"

    @staticmethod
    def base_command():
        return ["cat"]

    def inputs(self):
        return [
            ToolInput("files", Array(File()))
        ]

    def outputs(self):
        return [
            ToolOutput("merged", Stdout())
        ]
