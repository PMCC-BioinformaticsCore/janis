import unittest
from typing import List

from janis import ToolOutput, ToolInput, String, CommandTool, Stdout, InputSelector, Array, File
import janis.translations.cwl as cwl


class TestTool(CommandTool):

    @staticmethod
    def tool(): return "TestTranslation-tool"

    @staticmethod
    def base_command(): return "echo"

    def inputs(self) -> List[ToolInput]: return [ToolInput("testtool", String())]

    def outputs(self) -> List[ToolOutput]: return [ToolOutput("std", Stdout())]

    def friendly_name(self) -> str: return "Tool for testing translation"

    @staticmethod
    def docker(): return "ubuntu:latest"


class TestCwl(unittest.TestCase):

    def test_str_tool(self):
        t = TestTool()
        self.assertEqual(t.translate("cwl"), cwl_testtool)

    def test_input_selector_base(self):
        input_sel = InputSelector("random")
        self.assertEqual("$(inputs.random)", cwl.translate_input_selector(input_sel))

    def test_input_selector_prefix(self):
        input_sel = InputSelector("random", prefix="&& ")
        self.assertEqual("&& $(inputs.random)", cwl.translate_input_selector(input_sel))

    def test_base_input_selector(self):
        input_sel = InputSelector("random", suffix=".cwl")
        self.assertEqual("$(inputs.random).cwl", cwl.translate_input_selector(input_sel))


cwl_testtool = """\
baseCommand: echo
class: CommandLineTool
cwlVersion: v1.0
id: testtranslation-tool
inputs:
- id: testtool
  label: testtool
  type: string
label: testtranslation-tool
outputs:
- id: std
  label: std
  type: stdout
requirements:
  DockerRequirement:
    dockerPull: ubuntu:latest
  InlineJavascriptRequirement: {}
"""
