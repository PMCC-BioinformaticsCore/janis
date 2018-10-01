#
# Compile a java file
#

from typing import List, Dict

from examples.unix_commands.data_types import generic_file
from examples.unix_commands.data_types import class_file
from pipeline_definition.types.input_type import InputType
from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step


class Compile(Step):

  def provides(self) -> List[InputType]:
    return [class_file]

  def requires(self) -> List[InputType]:
    return [generic_file]

  def translate(self, mapped_inputs) -> Dict[str, str]:
    xlate = dict()

    xlate['run'] = '../tools/src/tools/compile.cwl'
    xlate['requirements'] = {'ResourceRequirement': {'coresMin': self.cores(), 'ramMin': self.ram()}}

    for mi in mapped_inputs:
      for candidate in mi.candidates.values():
        if mi.step_output_id == 'trimmed reads' and candidate['tag'] == self.tag():
          compile_step = candidate['step']

    inx = dict()

    inx['src'] = {'source': f'{compile_step}/tar_file'}
    inx['extractfile'] = 'hello.java'

    xlate['in'] = inx
    xlate['out'] = ['classfile']

    return {self.id(), xlate}

  def cores(self) -> int:
    return 2

  def ram(self) -> int:
    return 2*4000


class CompileFactory(StepFactory):
  @classmethod
  def type(cls) -> str:
    return 'compile'

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
  def build(cls, meta, debug=False) -> Compile:
    return Compile(meta, debug=debug)
