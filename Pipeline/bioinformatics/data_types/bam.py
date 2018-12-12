from Pipeline import File


class Bam(File):

    @staticmethod
    def name():
        return "BAM"

    @staticmethod
    def doc():
        return "A binary version of a SAM file, http://software.broadinstitute.org/software/igv/bam"


# import glob
#
# from pipeline_definition.types.input_type import InputFactory, InputType
# from pipeline_definition.types.input_type import Input
#
#
# bam_file_type = InputType('bam', label='a BAM file', description='A binary alignment map (BAM) file')
#
#
# class BamInput(Input):
#
#   def __init__(self, label: str, meta: dict):
#     super().__init__(label, meta)
#     self.path = None
#     self._files = []
#     self._resolved = False
#
#     if self.meta is not None:
#       self.path = self.meta().get("path")
#
#   def translate_for_workflow(self) -> dict:
#     raise Exception('translate_for_input not implemented for BAM files')
#
#   def translate_for_input(self) -> dict:
#     if self._resolved:
#       fd = [{'class': 'File', 'path': f} for f in self._files]
#     else:
#       fd = {'class': 'File', 'path': self.meta()['path']}
#
#     return {self.id(): fd}
#
#   def resolve(self):
#     self._resolved = Truex
#     self._files = glob.glob(self.meta()['path'])
#
#   def identify(self):
#     super().identify()
#
#   def datum_type(self) -> InputType:
#     return self.type()
#
#   def is_subtype_of(self, other) -> bool:
#     return False
#
#
# class BamFactory(InputFactory):
#   @classmethod
#   def type(cls) -> InputType:
#     return bam_file_type
#
#   @classmethod
#   def label(cls) -> str:
#     return bam_file_type.label()
#
#   @classmethod
#   def description(cls):
#     return bam_file_type.description()
#
#   @classmethod
#   def schema(cls):
#     return {
#       'schema': {
#         'path': {'type': 'string'},
#         'label': {'type': 'string'}
#       },
#       'nullable': True
#     }
#
#   @classmethod
#   def build(cls, label: str, meta: dict) -> BamInput:
#     return BamInput(label, meta)
