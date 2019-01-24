from typing import List

from Pipeline import ToolInput, File, ToolOutput, Array, String
from Pipeline.tool.expressiontool import ExpressionTool


class ScatterLines(ExpressionTool):
    @staticmethod
    def tool():
        return "parse-file"

    def expression(self):
        return """
${return { lines: inputs.file.contents.split("\n")
    .filter(function(q) { return q.length > 0; }) }}""".strip()

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("file", File())
        ]

    def outputs(self):
        return [ToolOutput("lines", Array(String()))]
