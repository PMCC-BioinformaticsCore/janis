#
# Untar a file
from typing import Dict

from examples.unix_commands.data_types.tar_file import tar_file
from pipeline_definition.types.step_type import Step
from pipeline_definition.types.step_type import StepFactory, StepInput, StepOutput


class Tar(Step):
    def requires(self) -> Dict[str, StepInput]:
        inp1 = self.get_input1()
        inp2 = self.get_input2()
        inp3 = self.get_input3()
        return {inp.tag: inp for inp in [inp1, inp2, inp3]}

    def provides(self) -> Dict[str, StepOutput]:
        outp = self.get_output()
        return {outp.tag: outp}

    def translate(self, mapped_inputs) -> dict:
        xlate = dict()

        xlate['run'] = '../tools/src/tools/tar.cwl'
        xlate['requirements'] = {'ResourceRequirement': {'coresMin': self.cores(), 'ramMin': self.ram()}}

        for mi in mapped_inputs:
            for candidate in mi.candidates.values():
                if mi.input_type == tar_file.type_name():  # and candidate['tag'] == self.tag():
                    tar_step = candidate['step']
                    tar_id = candidate['id']

        # inx = dict()
        #
        # inx['tarfile'] = {'source': f'{tar_step}/{tar_id}'}
        # inx['extractfile'] = 'hello.java'
        # xlate['in'] = inx
        # xlate['out'] = ['tar']

        return {self.id(): xlate}

    def cores(self) -> int:
        return 1

    def ram(self) -> int:
        return 1000

    @staticmethod
    def get_input1():
        return StepInput("input1", "File")

    @staticmethod
    def get_input2():
        return StepInput("input2", "File")

    @staticmethod
    def get_input3():
        return StepInput("tarName", "String")

    @staticmethod
    def get_output():
        return StepOutput("out", "File")


class TarFactory(StepFactory):
    @classmethod
    def type(cls) -> str:
        return 'tar'

    @classmethod
    def label(cls) -> str:
        return 'tar a file'

    @classmethod
    def description(cls) -> str:
        return 'tar 2 files and return the result.'

    @classmethod
    def schema(cls) -> dict:
        return {
            'schema': {
                'untar': {
                    'type': 'string'
                }
            },
            'nullable': True
        }

    @classmethod
    def build(cls, label: str, meta: dict) -> Tar:
        return Tar(label, meta)


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