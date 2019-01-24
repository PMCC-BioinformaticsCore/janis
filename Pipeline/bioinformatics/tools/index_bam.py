# from Workflow.step import StepFactory
# from Workflow.step import Step
#
#
# class IndexBamFactory(StepFactory):
#
#   @classmethod
#   def type(cls):
#     return 'index'
#
#   @classmethod
#   def label(cls):
#     return 'index a bam file'
#
#   @classmethod
#   def description(cls):
#     return 'index a bam file'
#
#   @classmethod
#   def schema(cls):
#     return {
#       'schema': {
#         'index': {}
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
#         Step.STR_ID: "indexedfile",
#         Step.STR_TYPE: "bamindex"
#       }
#     ]
#
#   def requires(self):
#     return [
#       {
#         Step.STR_ID: "sortedbamfile",
#         Step.STR_TYPE: "bam"
#       }
#     ]
#
#   def cores(self):
#     return 8
#
#   def ram(self):
#     return self.cores() * 4000
#
#   def translate(self, step_inputs):
#
#     xlate = dict()
#
#     xlate['run'] = '../tools/src/tools/samtools-index.cwl'
#     xlate['requirements'] = {'ResourceRequirement': {'coresMin': self.cores(), 'ramMin': self.ram()}}
#
#     for mi in step_inputs:
#       for candidate in mi.candidates.values():
#         if mi.step_output_id == 'sortedbamfile' and candidate['tag'] == self.tag():
#           align_step = candidate['step']
#
#     inx = dict()
#
#     inx['input'] = {'source': f'{align_step}/aligned-file'}
#     inx['threads'] = {'valueFrom': f'$( {self.cores()} )'}
#
#     xlate['in'] = inx
#     xlate['out'] = ['index']
#
#     return {self.id(): xlate}
