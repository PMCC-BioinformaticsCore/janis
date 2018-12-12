# from pipeline_definition.types.input_type import InputFactory, InputType
# from pipeline_definition.types.input_type import Input
#
# text_file_type = InputType('text', label='plain text')
#
#
# class TextFactory(InputFactory):
#   @classmethod
#   def type(cls) -> InputType:
#     return text_file_type
#
#   @classmethod
#   def label(cls):
#     return 'Text file'
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
#     return TextInput(input_dict, debug=debug)
#
#
# class TextInput(Input):
#   def translate_for_input(self):
#     pass
#
#   def __init__(self, input_dict, debug=False):
#     super().__init__(input_dict)
#     self.path = None
#     self._debug = debug
#
#     if self.meta is not None:
#       self.path = self.meta().get("path")
#
#   def identify(self):
#     super().identify()
#     print("Path:", self.path)
#
#   def datum_type(self):
#     return self.type()
#
#   def is_subtype_of(self, other):
#     return False
