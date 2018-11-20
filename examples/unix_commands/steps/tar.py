#
# Untar a file
from typing import List

from pipeline_definition.types.common_data_types import File, Number
from pipeline_definition.types.step import Tool, ToolInput, ToolOutput


class Tar(Tool):

    @staticmethod
    def tool():
        return "Tar"

    @staticmethod
    def supported_translations() -> List[str]:
        return ["cwl"]

    def inputs(self) -> List[ToolInput]:
        return [ToolInput("input1", File()), ToolInput("input2", File()), ToolInput("tarName", Number())]

    def outputs(self) -> List[ToolOutput]:
        return [ToolOutput("out", File())]



# class Tare(Step):
#     def requires(self) -> Dict[str, ToolInput]:
#         inp1 = self.get_input1()
#         inp2 = self.get_input2()
#         inp3 = self.get_input3()
#         return {inp.tag: inp for inp in [inp1, inp2, inp3]}
#
#     def provides(self) -> Dict[str, ToolOutput]:
#         outp = self.get_output()
#         return {outp.tag: outp}
#
#     def translate(self, mapped_inputs) -> dict:
#         xlate = dict()
#
#         xlate['run'] = '../tools/src/tools/tar.cwl'
#         xlate['requirements'] = {'ResourceRequirement': {'coresMin': self.cores(), 'ramMin': self.ram()}}
#
#         for mi in mapped_inputs:
#             for candidate in mi.candidates.values():
#                 if mi.input_type == tar_file.type_name():  # and candidate['tag'] == self.tag():
#                     tar_step = candidate['step']
#                     tar_id = candidate['id']
#
#         # inx = dict()
#         #
#         # inx['tarfile'] = {'source': f'{tar_step}/{tar_id}'}
#         # inx['extractfile'] = 'hello.java'
#         # xlate['in'] = inx
#         # xlate['out'] = ['tar']
#
#         return {self.id(): xlate}
#
#     def cores(self) -> int:
#         return 1
#
#     def ram(self) -> int:
#         return 1000
#
#     @staticmethod
#     def get_input1():
#         return ToolInput("input1", File)
#
#     @staticmethod
#     def get_input2():
#         return ToolInput("input2", File)
#
#     @staticmethod
#     def get_input3():
#         return ToolInput("tarName", String)
#
#     @staticmethod
#     def get_output():
#         return ToolOutput("out", TarFile)
#
#
# class TarFactory(StepFactory):
#     @classmethod
#     def type(cls) -> str:
#         return 'tar'
#
#     @classmethod
#     def label(cls) -> str:
#         return 'tar a file'
#
#     @classmethod
#     def description(cls) -> str:
#         return 'tar 2 files and return the result.'
#
#     @classmethod
#     def schema(cls) -> dict:
#         return {
#             'schema': {
#                 'untar': {
#                     'type': 'string'
#                 }
#             },
#             'nullable': True
#         }
#
#     @classmethod
#     def build(cls, label: str, meta: dict) -> Tar:
#         return Tar(label, meta)


"""

EXAMPLE XML

        <Step>
            <title>Tar</title>
            <author>Michael Franklin</author>
            <cores>1</cores>
            <ram>1024</ram>
            <label>{{label}}</label>
            <supportedTypes>
                <type>cwl</type>
                <type>wdl</type>
            </supportedTypes>
        
            <inputs>
                <input>
                    <tag>input2</tag>
                    <type>File</type>
                </input>
                <input>
                    <tag>input2</tag>
                    <type>File</type>
                </input>
            </inputs>
            <outputs>
                <output>
                    <tag>output</tag>
                    <type>File</type>
                </output>
            </outputs>
        </Step>

==============================

        class: CommandLineTool
        cwlVersion: v1.0
        id: tar
        baseCommand:
          - tar
          - cvf
        inputs:
          - id: tarName
            type: string?
            inputBinding:
              position: 0
          - id: input1
            type:
              - File
            inputBinding:
              position: 1
          - id: input2
            type:
              - File
            inputBinding:
              position: 2
        outputs:
          - id: output
            type: File?
            outputBinding:
              glob: '*'
        label: tar

  ========================================
  
        task SingleTar {
            String tarName
            File input1
            File input2
        
            command {
                tar cvf ${tarName} ${input1} ${input2}
            }
        
            output {
                File out = glob("*")[0]
            }
        }
"""