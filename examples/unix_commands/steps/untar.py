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
  def type(cls) -> InputType:
    pass

  @classmethod
  def label(cls) -> str:
    pass

  @classmethod
  def description(cls) -> str:
    pass

  @classmethod
  def schema(cls) -> dict:
    pass

  @classmethod
  def build(cls, meta, debug=False) -> TarFile:
    pass


class Untar(Step):
  def provides(self) -> List[InputType]:
    return [tar_file]

  def requires(self) -> List[InputType]:
    return [generic_file]

  def translate(self, mapped_inputs) -> dict:
    xlate = dict()

    xlate['run'] = '../tools/src/tools/tar-param.cwl'
    xlate['requirements'] = {'ResourceRequirement': {'coresMin': self.cores(), 'ramMin': self.ram()}}

    for mi in step_inputs:
      for candidate in mi.candidates.values():
        if mi.step_output_id == 'trimmed reads' and candidate['tag'] == self.tag():
          tar_step = candidate['step']
          read_id = candidate['id']
        if mi.step_output_id == 'reference':
          reference_step = candidate['step']
          reference_id = candidate['id']

    return xlate

  def cores(self) -> int:
    return 1

  def ram(self) -> int:
    return 1000
