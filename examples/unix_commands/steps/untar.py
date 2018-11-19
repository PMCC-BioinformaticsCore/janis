#
# Untar a file
#

from typing import Dict, Any

from examples.unix_commands.data_types.tar_file import tar_file
from pipeline_definition.types.step import Step
from pipeline_definition.types.step import StepFactory, ToolInput, ToolOutput
from pipeline_definition.types.types import TarFile, File


class Untar(Step):
    def requires(self) -> Dict[str, ToolInput]:
        inp = self.get_input1()
        return {inp.tag: inp}

    def provides(self) -> Dict[str, ToolOutput]:
        outp = self.get_output()
        return {outp.tag: outp}

    def translate(self, mapped_inputs) -> Dict[str, Any]:
        xlate: Dict[str, Any] = {
            'run': '../tools/src/tools/tar-param.cwl',
            'requirements': {
                'ResourceRequirement': {
                    'coresMin': self.cores(),
                    'ramMin': self.ram()
                }
            }
        }

        for mi in mapped_inputs:
            for candidate in mi.candidates.values():
                if mi.input_type == tar_file.type_name():  # and candidate['tag'] == self.tag():
                    tar_step = candidate['step']
                    tar_id = candidate['id']

        inx: Dict[str, Any] = {
            'tarfile': {
                'source': f'{tar_step}/{tar_id}'
            },
            'extractfile': 'hello.java'
        }

        xlate['in'] = inx
        xlate['out'] = ['untar']

        return {self.id(): xlate}

    def cores(self) -> int:
        return 1

    def ram(self) -> int:
        return 1000

    @staticmethod
    def get_input1():
        return ToolInput("input", TarFile)

    @staticmethod
    def get_output():
        return ToolOutput("out", File)


class UntarFactory(StepFactory):
    @classmethod
    def type(cls) -> str:
        return 'untar'

    @classmethod
    def label(cls) -> str:
        return 'untar a file'

    @classmethod
    def description(cls) -> str:
        return 'untar an archive and extract one or more files. Directories probably won\'t work.'

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
    def build(cls, label: str, meta: dict) -> Untar:
        return Untar(label, meta)


"""

EXAMPLE XML

<Step>
    <title>Untar</title>
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
            <tag>input</tag>
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

+ provide CWL template

#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [tar, xf]
inputs:
  input:
    type: File
    inputBinding:
      position: 1
  output:
    type: string
    inputBinding:
      position: 2

outputs:
  example_out:
    type: File
    outputBinding:
      glob: $(inputs.outputfile)

+ wdl template

task {
    File input
    File output
    command {
        tar -xf ${input}
    }
}

"""