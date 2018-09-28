from examples.bio_informatics.data_types.bam import *


class SortedBamFactory(BamFactory):
  @classmethod
  def type(cls):
    return 'sortedbam'

  @classmethod
  def build(cls, input_dict, debug=False):
    return SortedBamInput(input_dict, debug=debug)


class SortedBamInput(BamInput):

  def is_subtype_of(self, other):
    return isinstance(self, other)

