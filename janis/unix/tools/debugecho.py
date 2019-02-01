from typing import List

import janis.hints as hints
from janis import ToolOutput, ToolInput, String, Stdout, File
from janis.unix.tools.unixtool import UnixTool


class DebugEchoInputs(UnixTool):
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

    @staticmethod
    def hint_map():
        return {
            hints.CaptureType.KEY: {
                hints.CaptureType.TARGETED: {
                    "ramMin": 2000,
                    "coresMin": 3
                }
            }
        }
