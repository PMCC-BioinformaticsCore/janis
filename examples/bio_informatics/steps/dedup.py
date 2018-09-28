from typing import List

from pipeline_definition.types.input_type import InputType
from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step


class DedupFactory(StepFactory):
  @classmethod
  def type(cls):
    return 'dedup'

  @classmethod
  def label(cls):
    return 'dedup'

  @classmethod
  def description(cls):
    return cls.label()

  @classmethod
  def schema(cls):
    return {
      'schema': {
        'dedup': {
          'type': 'string'
        }
      },
      'nullable': True
    }

  @classmethod
  def build(cls, meta, debug=False):
    return DedupStep(meta, debug=debug)


class DedupStep(Step):

  def cores(self):
    return 2

  def ram(self):
    return 8000

  def translate(self, step_inputs):
    return {
        'command': 'dedup',
        'inputs': step_inputs
      }

  def provides(self) -> List[InputType]:
    return None
    # return [
    #   {
    #     Step.STR_ID: "bamfile",
    #     Step.STR_TYPE: "bam"
    #   }
    # ]

  def requires(self) -> List[InputType]:
    return None
    # return [
    #   {
    #     Step.STR_ID: "bamfile",
    #     Step.STR_TYPE: "bam"
    #   }
    # ]
