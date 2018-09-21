from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step


class CallFactory(StepFactory):
  @classmethod
  def type(cls):
    return 'call'

  @classmethod
  def label(cls):
    return 'call'

  @classmethod
  def description(cls):
    return cls.label()

  @classmethod
  def describe(cls):
    return {
      'schema': {
        'caller': {
          'type': 'string'
        }
      },
      'nullable': True
    }

  @classmethod
  def build(cls, meta, debug=False):
    step = CallStep(meta, debug=debug)
    return step


class CallStep(Step):

  def cores(self):
    return 1

  def ram(self):
    return 8000

  def translate(self, step_inputs):
    return {
        'command': 'vcf',
        'inputs': step_inputs
      }

  def provides(self):
    return [
      {
        Step.STR_ID: "read",
        Step.STR_TYPE: "VCF"
      }
    ]

  def requires(self):
    return [
      {
        Step.STR_ID: "alignedbamfile",
        Step.STR_TYPE: "bam"
      },
      {
        Step.STR_ID: "reference",
        Step.STR_TYPE: "reference"
      }
    ]
