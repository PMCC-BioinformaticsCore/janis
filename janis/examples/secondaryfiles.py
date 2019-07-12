from typing import List

from janis import ToolOutput, ToolInput, Workflow, Step, Input, Output, File, Array
from janis.types import InputSelector
from janis.unix.tools.unixtool import UnixTool


class DataTypeWithSecondary(File):
    def __init__(self, optional=False):
        super().__init__(optional, extension=".ext")

    @staticmethod
    def name():
        return "DataTypeWithSecondary"

    @staticmethod
    def secondary_files():
        return ["^.txt"]


class ToolThatAcceptsAndReturnsSecondary(UnixTool):
    @staticmethod
    def tool():
        return "ToolThatAcceptsAndReturnsSecondary"

    def friendly_name(self) -> str:
        return "Example tool only"

    @staticmethod
    def base_command():
        return "echo"  # non functional tool

    def inputs(self) -> List[ToolInput]:
        return [ToolInput("inp", DataTypeWithSecondary())]

    def outputs(self) -> List[ToolOutput]:
        return [ToolOutput("out", DataTypeWithSecondary(), glob=InputSelector("inp"))]


if __name__ == "__main__":
    w = Workflow("test_workflow")

    inp = Input("inp", DataTypeWithSecondary())
    stp = Step("stp", ToolThatAcceptsAndReturnsSecondary())
    w.add_pipe(inp, stp, Output("outp"))

    w.translate("wdl")

    w2 = Workflow("scattered_test_workflow")
    inp2 = Input("inp2", Array(DataTypeWithSecondary()), value=["path/to/file.ext"])
    stp2 = Step("stp2", ToolThatAcceptsAndReturnsSecondary())
    w2.add_pipe(inp2, stp2, Output("outp2"))
    w2.translate("wdl")
