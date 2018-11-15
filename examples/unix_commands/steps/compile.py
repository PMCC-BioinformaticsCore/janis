#
# Compile a java file
#

from typing import Dict

from examples.unix_commands.data_types.generic_file import generic_file
from pipeline_definition.types.step_type import Step, StepInput, StepOutput
from pipeline_definition.types.step_type import StepFactory


class Compile(Step):

    input1 = "input"

    def requires(self) -> Dict[str, StepInput]:
        inp1 = self.get_input1()
        return { inp1.tag: inp1 }

    def provides(self) -> Dict[str, StepOutput]:
        out = self.get_output()
        return { inp.tag: inp for inp in [out] }

    def translate(self, mapped_inputs) -> Dict[str, Dict]:
        xlate = dict()

        xlate['run'] = '../tools/src/tools/compile.cwl'
        xlate['requirements'] = {'ResourceRequirement': {'coresMin': self.cores(), 'ramMin': self.ram()}}

        for mi in mapped_inputs:
            for candidate in mi.candidates.values():
                if mi.input_type == generic_file.type_name() and candidate['tag'] == self.tag():
                    compile_step = candidate['step']
                    compile_id = candidate['id']

        inx = dict()

        inx['src'] = {'source': f'{compile_step}/{compile_id}'}
        inx['extractfile'] = 'hello.java'

        xlate['in'] = inx
        xlate['out'] = ['classfile']

        return {self.id(): xlate}

    def cores(self) -> int:
        return 2

    def ram(self) -> int:
        return 2 * 4000

    def get_input1(self) -> StepInput:
        return StepInput(Compile.input1, "File")

    def get_output(self) -> StepOutput:
        return StepOutput("out", "File")




class CompileFactory(StepFactory):
    @classmethod
    def type(cls) -> str:
        return 'java-compiler'

    @classmethod
    def label(cls) -> str:
        return 'compile a java file'

    @classmethod
    def schema(cls) -> dict:
        return {
            'schema': {
                'compile': {
                    'type': 'string'
                }
            },
            'nullable': True
        }

    @classmethod
    def build(cls, label: str, meta: dict) -> Compile:
        return Compile(label, meta)
