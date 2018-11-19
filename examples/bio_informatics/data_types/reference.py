# import os
#
# from pipeline_definition.utils.errors import PipelineTranslatorException
# from pipeline_definition.types.input_type import InputFactory, InputType
# from pipeline_definition.types.input_type import Input
#
# reference_file_type = InputType('reference', label='a reference genome')
#
#
# class ReferenceFactory(InputFactory):
#   @classmethod
#   def type(cls) -> InputType:
#     return reference_file_type
#
#   @classmethod
#   def label(cls):
#     return 'A reference genome'
#
#   @classmethod
#   def description(cls):
#     return cls.label()
#
#   @classmethod
#   def schema(cls):
#     return {
#       'schema': {
#         'path': {'type': 'string'},
#         'label': {'type': 'string'}
#       },
#       'nullable': True
#
#     }
#
#   @classmethod
#   def build(cls, input_dict, debug=False):
#     return ReferenceInput(input_dict)
#
#
# class ReferenceInput(Input):
#   def translate_for_input(self):
#     return {self.id(): {'class': 'File', 'path': self.meta()['path']}}
#
#   def translate_for_workflow(self):
#     return {self.id() + '_reference': 'File'}
#
#   def resolve(self):
#     if not os.path.exists(self.path):
#       raise PipelineTranslatorException('The reference {self.path} does not exist')
#
#   def __init__(self, input_dict, debug=False):
#     super().__init__(input_dict)
#     self.path = None
#     self.__debug = debug
#
#     if self.meta is not None:
#       self.path = self.meta().get("path")
#
#   def identify(self):
#     super().identify()
#     if self.__debug:
#       print("Path:", self.path)
#
#   def datum_type(self):
#     return self.type()
#
#   def is_subtype_of(self, other):
#     return False
