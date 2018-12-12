# import glob
#
# from pipeline_definition.utils.errors import PipelineTranslatorException
# from pipeline_definition.types.input_type import InputFactory, InputType
# from pipeline_definition.types.input_type import Input
#
# trimmed_reads_type = InputType('TrimmedReads', label='A pair of sequence files that have been trimmed')
#
#
# class TrimmedReads(Input):
#   def translate_for_input(self):
#     raise Exception('Translation noy yet implemented')
#
#   def translate_for_workflow(self):
#     raise Exception('Translation noy yet implemented')
#
#   def resolve(self):
#     self.forward_paired_files = glob.glob(self.forward_paired_pattern)
#     self.backward_paired_files = glob.glob(self.backward_paired_pattern)
#     self.forward_unpaired_files = glob.glob(self.forward_unpaired_pattern)
#     self.backward_unpaired_files = glob.glob(self.backward_unpaired_pattern)
#
#     if len(self.forward_paired_files) != len(self.backward_paired_files):
#       raise PipelineTranslatorException('The number of forward reads is different to the number of backward reads.')
#
#     self._resolved = True
#
#   def __init__(self, input_dict, debug=False):
#     super().__init__(input_dict, debug)
#     self._resolved = False
#     self.forward_paired_pattern = None
#     self.backward_paired_pattern = None
#     self.forward_unpaired_pattern = None
#     self.backward_unpaired_pattern = None
#     self.forward_paired_files = []
#     self.backward_paired_files = []
#     self.forward_unpaired_files = []
#     self.backward_unpaired_files = []
#     self.__debug = debug
#
#     if self.meta is not None:
#       self.forward_pattern = self.meta().get("forward-pattern")
#       self.backward_pattern = self.meta().get("backward-pattern")
#
#   def datum_type(self):
#     return self.type()
#
#   def is_subtype_of(self, other):
#     return False
#
#
# class TrimmedReadFactory(InputFactory):
#   @classmethod
#   def type(cls) -> InputType:
#     return trimmed_reads_type
#
#   @classmethod
#   def schema(cls) -> dict:
#     return {
#       'schema': {
#         'label': {'type': 'string'},
#         'paired-forward-pattern': {'type': 'string', 'required': True},
#         'paired-backward-pattern': {'type': 'string'},
#         'unpaired-forward-pattern': {'type': 'string', 'required': True},
#         'unpaired-backward-pattern': {'type': 'string'}
#       },
#       'nullable': True
#     }
#
#   @classmethod
#   def build(cls, input_dict, debug=False) -> TrimmedReads:
#     return TrimmedReads(input_dict, debug)
