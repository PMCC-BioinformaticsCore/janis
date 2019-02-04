# from typing import List
#
# from janis import ToolInput, File, ToolOutput, Array, String
# from janis.tool.expressiontool import ExpressionTool
# from janis.unix.tools.unixtool import UnixTool
#
#
# class ScatterLines(ExpressionTool, UnixTool):
#
#     @staticmethod
#     def base_command():
#         return None
#
#     def friendly_name(self) -> str:
#         return "Scatter across lines"
#
#     @staticmethod
#     def tool():
#         return "parse-file"
#
#     def expression(self):
#         return """
# ${return { lines: inputs.file.contents.split("\n")
#     .filter(function(q) { return q.length > 0; }) }}""".strip()
#
#     def inputs(self) -> List[ToolInput]:
#         return [
#             ToolInput("file", File())
#         ]
#
#     def outputs(self):
#         return [ToolOutput("lines", Array(String()))]
