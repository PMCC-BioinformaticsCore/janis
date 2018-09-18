
#
# Test discovery of inputs and generation of the input definition required by CWL.
import unittest
import yaml
from pipeline_definition.pipeline_translator import PipelineTranslator

import examples.bio_informatics


_yml = """
inputs:
  my_reads:
    SequenceReadArchivePaired:
      forward-pattern: '../test-data/*_R1.fastq.gz'
      backward-pattern: '../test-data/*_R2.fastq.gz'
  my_ref:
    reference:
      path: ../test-data/hg38_no_alt.fa
  my_bam:
    bam:
      path: ../test-data/aligned.bam
"""

_expected_resolved = yaml.load("""
inputs:
  my_reads_backward:
  - class: File
    path: ../test-data/S1_R2.fastq.gz
  - class: File
    path: ../test-data/S2_R2.fastq.gz
  my_reads_forward:
  - class: File
    path: ../test-data/S1_R1.fastq.gz
  - class: File
    path: ../test-data/S2_R1.fastq.gz
  my_ref:
    class: File
    path: ../test-data/hg38_no_alt.fa
  my_bam:
    - class: File
      path: ../test-data/aligned.bam
""")

_expected_unresolved = yaml.load("""
inputs:
  my_reads_backward:
    class: File
    path: ../test-data/*_R2.fastq.gz
  my_reads_forward:
    class: File
    path: ../test-data/*_R1.fastq.gz
  my_ref:
    class: File
    path: ../test-data/hg38_no_alt.fa
  my_bam:
    class: File
    path: ../test-data/aligned.bam
""")


class InputFileTest(unittest.TestCase):

  def translate(self, resolve, expected):
    translator = PipelineTranslator(debug=False)
    translator.translate_string(_yml)

    translation = translator.input(resolve=resolve)
    got = yaml.load(translation)

    self.assertTrue(got == expected)

  def test_resolved(self):
    self.translate(True, _expected_resolved)

  def test_unresolved(self):
    self.translate(False, _expected_unresolved)


if __name__ == '__main__':
  unittest.main()
