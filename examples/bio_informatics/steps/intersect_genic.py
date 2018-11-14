from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step, StepInput, StepOutput

from typing import List, Dict

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
  def schema(cls):
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
    return IntersectStep(meta)


class IntersectStep(Step):

    def translate(self, mapped_inputs) -> Dict[str, str]:
        pass

    def cores(self) -> int:
        return 1

    def ram(self) -> int:
        return 512

    def requires(self) -> Dict[str, StepInput]:
        inp = self.get_input1()
        return { inp.tag: inp }

    def provides(self) -> Dict[str, StepOutput]:
        outp = self.get_output()
        return { outp.tag: outp }

    def get_input1(self) -> StepInput:
        return StepInput("read", "SequenceReadArchivePaired")

    def get_output(self) -> StepOutput:
        return StepOutput("reports", "Text")

    # def provides(self):
    #   return [
    #     {
    #       Step.STR_ID: "reports",
    #       Step.STR_TYPE: "Text"
    #     }
    #   ]
    #
    # def requires(self):
    #   return [
    #     {
    #       Step.STR_ID: "read",
    #       Step.STR_TYPE: "SequenceReadArchivePaired"
    #     }
    #   ]
