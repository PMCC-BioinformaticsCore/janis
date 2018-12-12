# from Workflow.step import StepFactory
# from Workflow.step import Step
#
#
# class JointCallFactory(StepFactory):
#   @classmethod
#   def type(cls):
#     return 'joint_call'
#
#   @classmethod
#   def label(cls):
#     return 'joint_call'
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
#           'allowed': ['mutect'],
#           'default': 'mutect'
#         },
#         'normal_tag': {
#           'type': 'string'
#         },
#         'tumour_tag': {
#           'type': 'string'
#         }
#       },
#       'nullable': True
#     }
#
#   @classmethod
#   def build(cls, meta, debug=False):
#     step = JointCallStep(meta, debug=debug)
#     return step
#
#
# class JointCallStep(Step):
#
#   def cores(self):
#     return 8
#
#   def ram(self):
#     return 16000
#
#   def translate(self, step_inputs):
#     return {
#         'command': 'joint-call',
#         'inputs': step_inputs
#       }
#
#   def provides(self):
#     return [
#       {
#         Step.STR_ID: 'joint calls',
#         Step.STR_TYPE: 'vcf'
#       }
#     ]
#
#   def requires(self):
#     return [
#       {
#         Step.STR_ID: 'normal_tag',
#         Step.STR_TYPE: 'sortedbam'
#       },
#       {
#         Step.STR_ID: 'tumour_tag',
#         Step.STR_TYPE: 'sortedbam'
#       },
#       {
#         Step.STR_ID: 'reference',
#         Step.STR_TYPE: 'reference'
#       }
#     ]
#
#   def translate(self, step_inputs):
#     xlate = dict()
#
#     xlate['run'] = '../tools/src/tools/platypus.cwl'
#     xlate['requirements'] = {'ResourceRequirement': {'coresMin': self.cores(), 'ramMin': self.ram()}}
#
#     for mi in step_inputs:
#       for candidate in mi.candidates.values():
#         if mi.step_output_id == 'tumour_tag' and candidate['tag'] == 'tumour':
#           tumour_step = candidate['step']
#         if mi.step_output_id == 'normal_tag' and candidate['tag'] == 'normal':
#           normal_step = candidate['step']
#         if mi.step_output_id == 'reference':
#           reference_step = candidate['step']
#           reference_id = candidate['id']
#
#     inx = dict()
#
#     inx['tumour'] = {'source': f'{tumour_step}/sorted'}
#     inx['normal'] = {'source': f'{normal_step}/sorted'}
#     inx['reference'] = {'source': f'{reference_step}/{reference_id}'}
#     inx['vcf'] = {
#       'source': f'{tumour_step}/sorted',
#       'valueFrom': f'${{ return self.nameroot + ".vcf"; }}'
#     }
#
#     xlate['in'] = inx
#     xlate['out'] = ['coverage', 'call_stats', 'mutations']
#
#     return {self.id(): xlate}
