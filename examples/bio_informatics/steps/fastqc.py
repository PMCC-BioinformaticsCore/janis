from pipeline_definition.types.step import StepFactory, Step, ToolInput, ToolOutput
from pipeline_definition.types.step import Step
from pipeline_definition.types.input_type import InputType

from typing import List, Dict


class FastQCFactory(StepFactory):
    @classmethod
    def type(cls):
        return 'fastqc'

    @classmethod
    def label(cls):
        return 'fastqc'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def schema(cls):
        return {
            'schema': {},
            'nullable': True
        }

    @classmethod
    def build(cls, meta, debug=False):
        return FastQCStep(meta, debug=debug)



class FastQCStep(Step):

    def translate(self):
        pass

    def input_labels(self) -> [str]:
        return [FastQCStep.input1]

    def requires(self) -> Dict[str, ToolInput]:
        inp = self.get_input1()
        return { inp.tag: inp }

    def provides(self) -> Dict[str, ToolOutput]:
        outp = self.get_output()
        return { outp.tag: outp }

    # def provides(self):
    #     return [
    #         {
    #             Step.STR_ID: "reports",
    #             Step.STR_TYPE: "Text"
    #         }
    #     ]
    #
    # def requires(self) -> List[StepInput]:
    #     """
    #     Only returns a SequenceReadArchivePaired object, with tag: input
    #     :return:
    #     """
    #     return [self.get_input1()]

    def get_input1(self) -> ToolInput:
        return ToolInput("input", "SequenceReadArchivePaired")

    def get_output(self) -> ToolOutput:
        return ToolOutput("reports", "Text")
