from typing import List

from janis.types import WildcardSelector, InputSelector

from janis import ToolOutput, ToolInput, Workflow, Step, Input, Output, File, Array
from janis.tool.commandtool import CommandTool


class DataTypeWithSecondary(File):
    @staticmethod
    def name():
        return "DataTypeWithSecondary"

    @staticmethod
    def secondary_files():
        return [".txt"]


class ToolThatAcceptsAndReturnsSecondary(CommandTool):
    @staticmethod
    def tool():
        return "ToolThatAcceptsAndReturnsSecondary"

    def friendly_name(self) -> str:
        return "TestTool: Do not use"

    @staticmethod
    def base_command():
        return "echo" # non functional tool

    def inputs(self) -> List[ToolInput]:
        return [ToolInput("input", DataTypeWithSecondary())]

    def outputs(self) -> List[ToolOutput]:
        return [ToolOutput("output", DataTypeWithSecondary(), glob=InputSelector("input"))]


if __name__ == "__main__":
    w = Workflow("test_workflow")

    inp = Input("inp", DataTypeWithSecondary())
    stp = Step("stp", ToolThatAcceptsAndReturnsSecondary())
    w.add_pipe(inp, stp, Output("outp"))

    w.translate("wdl")

    w2 = Workflow("scattered_test_workflow")
    inp2 = Input("inp2", Array(DataTypeWithSecondary()))
    stp2 = Step("stp2", ToolThatAcceptsAndReturnsSecondary())
    w2.add_pipe(inp2, stp2, Output("outp2"))
    w2.translate("wdl")

