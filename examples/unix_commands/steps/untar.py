#
# Untar a file
#

from typing import List, Dict

from examples.unix_commands.data_types.generic_file import generic_file
from examples.unix_commands.data_types.tar_file import tar_file, TarFile
from pipeline_definition.types.input_type import InputType
from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step


class Untar(Step):
  def provides(self) -> Dict[str, InputType]:
    return {'untar': generic_file}

  def requires(self) -> List[InputType]:
    return [tar_file]

  def translate(self, mapped_inputs) -> dict:
    xlate = dict()

    xlate['run'] = '../tools/src/tools/tar-param.cwl'
    xlate['requirements'] = {'ResourceRequirement': {'coresMin': self.cores(), 'ramMin': self.ram()}}

    for mi in mapped_inputs:
      for candidate in mi.candidates.values():
        if mi.input_type == tar_file.type_name():  # and candidate['tag'] == self.tag():
          tar_step = candidate['step']
          tar_id = candidate['id']

    inx = dict()

    inx['tarfile'] = {'source': f'{tar_step}/{tar_id}'}
    inx['extractfile'] = 'hello.java'

    xlate['in'] = inx
    xlate['out'] = ['untar']

    return {self.id(): xlate}

  def cores(self) -> int:
    return 1

  def ram(self) -> int:
    return 1000


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
  def build(cls, label: str, meta: dict, debug=False) -> Untar:
    return Untar(label, meta)
