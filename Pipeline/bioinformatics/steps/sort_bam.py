# from Workflow.step import StepFactory
# from Workflow.step import Step
#
#
# class SortBamFactory(StepFactory):
#
#   @classmethod
#   def type(cls):
#     return 'sort'
#
#   @classmethod
#   def label(cls):
#     return 'sort a bam file'
#
#   @classmethod
#   def description(cls):
#     return 'sort a bam file'
#
#   @classmethod
#   def schema(cls):
#     return {
#       'schema': {
#         'sort': {
#           'type': 'string'
#         }
#       },
#       'nullable': True
#     }
#
#   @classmethod
#   def build(cls, meta, debug=False):
#     return SortBam(meta, debug=debug)
#
#
# class SortBam(Step):
#   def provides(self):
#     return [
#       {
#         Step.STR_ID: "sortedbamfile",
#         Step.STR_TYPE: "sortedbam"
#       }
#     ]
#
#   def requires(self):
#     return [
#       {
#         Step.STR_ID: "bamfile",
#         Step.STR_TYPE: "bam"
#       }
#     ]
#
#
#   def cores(self):
#     return 8
#
#   def ram(self):
#     return self.cores()*2000
#
#   def translate(self, step_inputs):
#
#     xlate = dict()
#
#     xlate['run'] = '../tools/src/tools/samtools-sort.cwl'
#     xlate['requirements'] = {'ResourceRequirement': {'coresMin': self.cores(), 'ramMin': self.ram()}}
#
#     for mi in step_inputs:
#       for candidate in mi.candidates.values():
#         if mi.step_output_id == 'bamfile' and candidate['tag'] == self.tag():
#           align_step = candidate['step']
#
#     inx = dict()
#
#     inx['input'] = {'source': f'{align_step}/aligned-file'}
#     inx['output_name'] = {
#       'source': f'{align_step}/aligned-file',
#       'valueFrom': f'${{return self.nameroot + ".sorted.{self.tag()}.bam";}}'
#     }
#     inx['threads'] = {'valueFrom': f'$( {self.cores()} )'}
#
#     xlate['in'] = inx
#     xlate['out'] = ['sorted']
#
#     return {self.id(): xlate}
