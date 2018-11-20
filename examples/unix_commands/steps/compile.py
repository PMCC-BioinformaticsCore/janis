#
# Compile a java file
#

from typing import List

from pipeline_definition.types.common_data_types import File
from pipeline_definition.types.step import ToolInput, ToolOutput, Tool


class Compile(Tool):

    @staticmethod
    def tool():
        return "java-compiler"

    @staticmethod
    def supported_translations() -> List[str]:
        return ["cwl"]

    def inputs(self) -> List[ToolInput]:
        return [ToolInput("input", File())]

    def outputs(self) -> List[ToolOutput]:
        return [ToolOutput("out", File())]


# class Compilee(Step):
#
#     input1 = "input"
#
#     def requires(self) -> Dict[str, ToolInput]:
#         inp1 = self.get_input1()
#         return { inp1.tag: inp1 }
#
#     def provides(self) -> Dict[str, ToolOutput]:
#         out = self.get_output()
#         return { inp.tag: inp for inp in [out] }
#
#     def translate(self, mapped_inputs) -> Dict[str, Any]:
#         xlate: Dict[str, Any] = {
#             'run': '../tools/src/tools/compile.cwl',
#             'requirements': {
#                 'ResourceRequirement': {
#                     'coresMin': self.cores(),
#                     'ramMin': self.ram()
#                 }
#             }
#         }
#
#         for mi in mapped_inputs:
#             for candidate in mi.candidates.values():
#                 if mi.input_type == generic_file.type_name() and candidate['tag'] == self.tag():
#                     compile_step = candidate['step']
#                     compile_id = candidate['id']
#
#         inx: Dict[str, Any] = {
#             'src': {
#                 'source': f'{compile_step}/{compile_id}'
#             },
#             'extractfile': 'hello.java'
#         }
#
#         xlate['in'] = inx
#         xlate['out'] = ['classfile']
#
#         return {self.id(): xlate}
#
#     def cores(self) -> int:
#         return 2
#
#     def ram(self) -> int:
#         return 2 * 4000
#
#     def get_input1(self) -> ToolInput:
#         return ToolInput(Compile.input1, "File")
#
#     def get_output(self) -> ToolOutput:
#         return ToolOutput("out", "File")
#
#
# class CompileFactory(StepFactory):
#     @classmethod
#     def type(cls) -> str:
#         return 'java-compiler'
#
#     @classmethod
#     def label(cls) -> str:
#         return 'compile a java file'
#
#     @classmethod
#     def schema(cls) -> dict:
#         return {
#             'schema': {
#                 'compile': {
#                     'type': 'string'
#                 }
#             },
#             'nullable': True
#         }
#
#     @classmethod
#     def build(cls, label: str, meta: dict) -> Compile:
#         return Compile(label, meta)
