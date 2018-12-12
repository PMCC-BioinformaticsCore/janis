# from Workflow.step import StepFactory
# from Workflow.step import Step
#
#
# class CallFactory(StepFactory):
#   @classmethod
#   def type(cls):
#     return 'call'
#
#   @classmethod
#   def label(cls):
#     return 'call'
#
#   @classmethod
#   def description(cls):
#     return cls.label()
#
#   @classmethod
#   def schema(cls):
#     return {
#       'schema': {
#         'caller': {
#           'type': 'string',
#           'allowed': ['platypus', 'vep'],
#           'default': 'platypus'
#         }
#       },
#       'nullable': True
#     }
#
#   @classmethod
#   def build(cls, meta, debug=False):
#     step = CallStep(meta, debug=debug)
#     return step
#
#
# class CallStep(Step):
#
#   def __init__(self, meta, debug=False):
#     super().__init__(meta, debug=debug)
#     meta = self.meta()
#     if meta is not None and meta['caller'] != 'platypus':
#       raise Exception('Sorry, only platypus is supported at the moment')
#
#   def cores(self):
#     return 2
#
#   def ram(self):
#     return 16000
#
#   def translate(self, step_inputs):
#     xlate = dict()
#
#     xlate['run'] = '../tools/src/tools/platypus.cwl'
#     xlate['requirements'] = {'ResourceRequirement': {'coresMin': self.cores(), 'ramMin': self.ram()}}
#
#     # inx = dict()
#     #
#     # inx['bamFiles'] = {
#     #   'source':
#     # }
#     return {self.id(): xlate}
#
#   def provides(self):
#     return [
#       {
#         Step.STR_ID: "calls",
#         Step.STR_TYPE: "vcf"
#       }
#     ]
#
#   def requires(self):
#     return [
#       {
#         Step.STR_ID: "indexedbam",
#         Step.STR_TYPE: "bam"
#       },
#       {
#         Step.STR_ID: "reference",
#         Step.STR_TYPE: "reference"
#       }
#     ]
