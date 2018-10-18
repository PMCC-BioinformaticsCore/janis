#
# Compile a java file
#

from typing import List, Dict

from examples.unix_commands.data_types.generic_file import generic_file
from examples.unix_commands.data_types.class_file import class_file
from pipeline_definition.types.input_type import InputType
from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step


class Compile(Step):

  def provides(self) -> List[InputType]:
    return [class_file]

  def requires(self) -> List[InputType]:
    return [generic_file]

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
