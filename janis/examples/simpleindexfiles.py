from typing import List
from janis import ToolOutput, ToolInput, Workflow, Step, Input, Output, File
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
        return [ToolOutput("output", DataTypeWithSecondary())]


if __name__ == "__main__":
    w = Workflow("test_workflow")

    inp = Input("inp", DataTypeWithSecondary())
    stp = Step("stp", ToolThatAcceptsAndReturnsSecondary())
    w.add_pipe(inp, stp, Output("outp"))

    w.dump_translation("wdl")