import unittest

from pipeline_definition.pipeline_translator import PipelineTranslator
import json
import examples.bio_informatics


_yml = """
inputs:
  fastq:
    SequenceReadArchivePaired:
      forward-pattern: '*_R1.fastq.gz'
      backward-pattern: '*_R2.fastq.gz'
  ref:
    bam:
      path: 'path/to/reference'

steps:
  - step1:
      trim:
        trimmer : 'trimmomatic'
  - step2:
      input_scope: ['ref']
      align:
        aligner: 'bwa'
  - step3:
      dedup:
  - step4:
      input_scope: ['step2', 'step3']
      call:
"""

_expected = json.loads("""""")


class LinearPipeline(unittest.TestCase):

  def test_graph(self):
    translator = PipelineTranslator(debug=False)
    translator.translate_string(_yml)
    translation = translator.pipeline()
    tr_json = json.loads(translation)
    print('-'*80)
    print(translation)
    print('-'*80)
    self.assertTrue(translation == _expected)


if __name__ == '__main__':
    unittest.main()
