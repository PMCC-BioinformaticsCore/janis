from typing import List

from Pipeline import CommandTool, ToolOutput, ToolInput, String, Stdout, File


class DebugEchoInputs(CommandTool):
    @staticmethod
    def tool():
        return "debuginputs"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self) -> List[ToolInput]:
        return [
            # add whatever you want here
            ToolInput("stringInput", String(optional=True)),
            ToolInput("fileInput", File(optional=True))
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("output", Stdout())
        ]

    def friendly_name(self) -> str:
        return "Debug inputs (log)"
