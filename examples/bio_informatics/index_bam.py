from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step


class IndexBamFactory(StepFactory):

  @classmethod
  def type(cls):
    return 'index'

  @classmethod
  def label(cls):
    return 'index a bam file'

  @classmethod
  def description(cls):
    return 'index a bam file'

  @classmethod
  def describe(cls):
    return {
      'schema': {
        'index': {}
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
        Step.STR_ID: "indexedfile",
        Step.STR_TYPE: "bamindex"
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
    return self.cores() * 8000
