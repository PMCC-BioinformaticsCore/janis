from typing import List

from janis_core import (
    ToolOutput,
    ToolInput,
    WorkflowBuilder,
    File,
    Array,
    ScatterDescription,
    ScatterMethods,
)
from janis_core.types import InputSelector
from janis_core.utils.scatter import ScatterMethod
from janis_unix.tools.unixtool import UnixTool
from janis_bioinformatics.data_types import FastaBwa, BamBai


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


class ToolTypeThatAcceptsMultipleBioinfTypes(UnixTool):
    def tool(self):
        return "TESTTOOL_BIOINF"

    def inputs(self):
        return [ToolInput("bam", BamBai), ToolInput("reference", FastaBwa)]

    def outputs(self):
        return [
            ToolOutput("out_bam", BamBai, glob=InputSelector("bam")),
            ToolOutput("out_reference", FastaBwa, glob=InputSelector("reference")),
        ]

    def base_command(self):
        return "echo"


if __name__ == "__main__":
    w = WorkflowBuilder("test_workflow")

    # EXAMPLE 1

    w.input("inp", DataTypeWithSecondary)
    w.step("stp", ToolThatAcceptsAndReturnsSecondary(inp=w.inp))
    w.output("out", source=w.stp)
    w.translate("wdl")

    # EXAMPLE 2

    w2 = WorkflowBuilder("scattered_test_workflow")
    w2.input("inp", Array(DataTypeWithSecondary), value=["path/to/file.ext"])
    w2.step("stp", ToolThatAcceptsAndReturnsSecondary(inp=w2.inp), scatter="inp")
    w2.output("out", source=w2.stp)
    w2.translate("wdl")

    # EXAMPLE 3

    w3 = WorkflowBuilder("scattered_bioinf_complex")
    w3.input("my_bams", Array(BamBai))
    w3.input("my_references", Array(FastaBwa))

    w3.step(
        "my_step",
        ToolTypeThatAcceptsMultipleBioinfTypes(
            bam=w3.my_bams, reference=w3.my_references
        ),
        scatter=ScatterDescription(["bam", "reference"], ScatterMethods.cross),
    )

    w3.output("out_bam", source=w3.my_step.out_bam)
    w3.output("out_reference", source=w3.my_step.out_reference)

    w3.translate("wdl")
