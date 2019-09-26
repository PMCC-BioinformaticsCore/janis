from typing import List

from janis_core import ToolOutput, ToolInput, WorkflowBuilder, File, Array
from janis_core.types import InputSelector
from janis_unix.tools.unixtool import UnixTool


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
    w = WorkflowBuilder("test_workflow")

    w.input("inp", DataTypeWithSecondary)
    w.step("stp", ToolThatAcceptsAndReturnsSecondary(inp=w.inp))
    w.output("out", source=w.stp)
    w.translate("wdl")

    w2 = WorkflowBuilder("scattered_test_workflow")
    w2.input("inp", Array(DataTypeWithSecondary), default=["path/to/file.ext"])
    w2.step("stp", ToolThatAcceptsAndReturnsSecondary(inp=w2.inp), scatter="inp")
    w2.output("out", source=w2.stp)
    w2.translate("wdl")
