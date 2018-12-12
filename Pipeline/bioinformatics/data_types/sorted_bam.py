# from examples.bio_informatics.data_types.bam import *
#
#
# sorted_bam_file_type = InputType('sorted_bam',
#                                  label='a sorted BAM file',
#                                  description='A binary alignment map (BAM) output from samtools sort  ')
#
#
# class SortedBamFactory(BamFactory):
#   @classmethod
#   def type(cls) -> InputType:
#     return sorted_bam_file_type
#
#   @classmethod
#   def build(cls, input_dict, debug=False):
#     return SortedBamInput(input_dict, debug=debug)
#
#
# class SortedBamInput(BamInput):
#
#   def is_subtype_of(self, other):
#     return isinstance(self, other)
#
