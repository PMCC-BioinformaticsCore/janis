import unittest
from typing import List

from janis import ToolOutput, ToolInput, String, CommandTool, Stdout


class TestTool(CommandTool):

    @staticmethod
    def tool(): return "TestTranslation-tool"

    @staticmethod
    def base_command(): return "echo"

    def inputs(self) -> List[ToolInput]: return [ToolInput("testtool", String())]

    def outputs(self) -> List[ToolOutput]: return [ToolOutput("std", Stdout())]

    def friendly_name(self) -> str: return "Tool for testing translation"


class TestTranslation(unittest.TestCase):

    def test_str_tool(self):
        t = TestTool()
        self.assertEqual(t.translate("cwl"), cwl_testtool)


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
