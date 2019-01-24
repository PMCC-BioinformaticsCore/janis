# from Workflow.step import StepFactory
# from Workflow.step import Step, ToolInput, ToolOutput
#
# from typing import Dict
#
# class DedupFactory(StepFactory):
#   @classmethod
#   def type(cls):
#     return 'dedup'
#
#   @classmethod
#   def label(cls):
#     return 'dedup'
#
#   @classmethod
#   def description(cls):
#     return cls.label()
#
#   @classmethod
#   def schema(cls):
#     return {
#       'schema': {
#         'dedup': {
#           'type': 'string'
#         }
#       },
#       'nullable': True
#     }
#
#   @classmethod
#   def build(cls, meta, debug=False):
#     return DedupStep(meta)
#
#
# class DedupStep(Step):
#
#   def cores(self):
#     return 2
#
#   def ram(self):
#     return 8000
#
#   def translate(self, step_inputs):
#     return {
#         'command': 'dedup',
#         'inputs': step_inputs
#       }
#
#   # def provides(self) -> List[InputType]:
#   #   return None
#   #   # return [
#   #   #   {
#   #   #     Step.STR_ID: "bamfile",
#   #   #     Step.STR_TYPE: "bam"
#   #   #   }
#   #   # ]
#   #
#   # def requires(self) -> List[InputType]:
#   #   return None
#   #   # return [
#   #   #   {
#   #   #     Step.STR_ID: "bamfile",
#   #   #     Step.STR_TYPE: "bam"
#   #   #   }
#   #   # ]
#
#   def provides(self) -> Dict[str, ToolOutput]:
#     # return [trimmed_reads_type]
#     outp = self.get_output()
#     return {outp.tag: outp}
#
#   def requires(self) -> Dict[str, ToolInput]:
#     # return [paired_reads_type]
#     inp = self.get_input()
#     return {inp.tag: inp}
#
#   def get_input(self) -> ToolInput:
#     return ToolInput("bamfile", "bam")
#
#   def get_output(self) -> ToolOutput:
#     return ToolOutput("bamfile", "bam")
