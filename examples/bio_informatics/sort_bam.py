from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step

class SortBamFactory(StepFactory):

  @classmethod
  def type(cls):
    return 'sort'

  @classmethod
  def label(cls):
    return 'sort a bam file'

  @classmethod
  def description(cls):
    return 'sort a bam file'

  @classmethod
  def describe(cls):
    return {
      'schema': {
        'sort': {
          'type': 'string'
        }
      },
      'nullable': True
    }

  @classmethod
  def build(cls, meta, debug=False):
    return SortBam(meta, debug=debug)


class SortBam(Step):
  def provides(self):
    return [
      {
        Step.STR_ID: "sortedfile",
        Step.STR_TYPE: "sortedbam"
      }
    ]

  def requires(self):
    return [
      {
        Step.STR_ID: "bamfile",
        Step.STR_TYPE: "bam"
      }
    ]

  def translate(self, mapped_inputs):
    pass

  def cores(self):
    return 2

  def ram(self):
    return self.cores()*8000
