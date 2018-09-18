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
  def describe(cls):
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

  def translate(self):
    pass

  def provides(self):
    return [
      {
        Step.STR_ID: "bamfile",
        Step.STR_TYPE: "bam"
      }
    ]

  def requires(self):
    return [
      {
        Step.STR_ID: "bamfile",
        Step.STR_TYPE: "bam"
      }
    ]
