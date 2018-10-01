#
# Untar a file
#

from typing import List

from examples.unix_commands.data_types import generic_file
from examples.unix_commands.data_types.tar_file import tar_file, TarFile
from pipeline_definition.types.input_type import InputType
from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step


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
  def build(cls, meta, debug=False) -> TarFile:
    return Untar(meta, debug=debug)


class Untar(Step):
  def provides(self) -> List[InputType]:
    return [tar_file]

  def requires(self) -> List[InputType]:
    return [generic_file]

  def translate(self, mapped_inputs) -> dict:
    xlate = dict()

    xlate['run'] = '../tools/src/tools/tar-param.cwl'
    xlate['requirements'] = {'ResourceRequirement': {'coresMin': self.cores(), 'ramMin': self.ram()}}

    for mi in mapped_inputs:
      for candidate in mi.candidates.values():
        if mi.step_output_id == 'trimmed reads' and candidate['tag'] == self.tag():
          tar_step = candidate['step']

    inx = dict()

    inx['tarfile'] = {'source': f'{tar_step}/tar_file'}
    inx['extractfile'] = 'hello.java'

    xlate['in'] = inx
    xlate['out'] = ['untar']

    return {self.id(), xlate}

  def cores(self) -> int:
    return 1

  def ram(self) -> int:
    return 1000
