import unittest

from pipeline_definition.pipeline_translator import PipelineTranslator
import examples.bio_informatics

_yml = """
inputs:
  fastq:
    SequenceReadArchivePaired:
      forward-pattern: '*_R1.fastq.gz'
      backward-pattern: '*_R2.fastq.gz'
  ref:
    REFERENCE:
      path: 'path/to/reference'

steps:
  - qc:
      fastqc:
  - remove-adaptors:
      trim:
  - align-to-human:
      align:
  - dedup:
      dedup:
  - intersect-genic:
      input_scope: [dedup]
      bedtools-intersect:
        split: true
  - intersect-nongenic:
      input_scope: [dedup]
      bedtools-intersect:
        reportNoOverlaps: true
"""


class BranchedPipeline(unittest.TestCase):

  def test_graph(self):
    translator = PipelineTranslator(debug=True)
    translation = translator.translate_string(_yml)
    # print(translation)
    self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()


