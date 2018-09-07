from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step


class IntersectFactory(StepFactory):
  @classmethod
  def type(cls):
    return 'bedtools-intersect'

  @classmethod
  def label(cls):
    return 'bedtools-intersect'

  @classmethod
  def description(cls):
    return cls.label()

  @classmethod
  def describe(cls):
    return {
      'schema': {
        "split": {
          "type": "boolean",
          "default": False
        },
        "reportNoOverlaps": {
          "type": "boolean",
          "default": False
        }
      },
      'nullable': True
    }

  @classmethod
  def build(cls, meta, debug=False):
    return IntersectStep(meta, debug=debug)


class IntersectStep(Step):

  def provides(self):
    return [
      {
        Step.STR_ID: "reports",
        Step.STR_TYPE: "Text"
      }
    ]

  def requires(self):
    return [
      {
        Step.STR_ID: "read",
        Step.STR_TYPE: "SequenceReadArchivePaired"
      }
    ]
