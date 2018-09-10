from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step


class JointCallFactory(StepFactory):
  @classmethod
  def type(cls):
    return 'joint_call'

  @classmethod
  def label(cls):
    return 'joint_call'

  @classmethod
  def description(cls):
    return cls.label()

  @classmethod
  def describe(cls):
    return {
      'schema': {
        'caller': {
          'type': 'string'
        },
        'normal_tag': {
          'type': 'string'
        },
        'tumour_tag': {
          'type': 'string'
        }
      },
      'nullable': True
    }

  @classmethod
  def build(cls, meta, debug=False):
    step = JointCallStep(meta, debug=debug)
    return step


class JointCallStep(Step):

  def provides(self):
    return [
      {
        Step.STR_ID: 'out1',
        Step.STR_TYPE: 'VCF'
      }
    ]

  def requires(self):
    return [
      {
        Step.STR_ID: 'normal_tag',
        Step.STR_TYPE: 'BAM'
      },
      {
        Step.STR_ID: 'tumour_tag',
        Step.STR_TYPE: 'BAM'
      },
      {
        Step.STR_ID: 'references',
        Step.STR_TYPE: 'REFERENCE'
      }
    ]
