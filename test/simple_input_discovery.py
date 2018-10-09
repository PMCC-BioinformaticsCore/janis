
#
# Test discovery of inputs in the file system and generation of the input definition required by CWL.
#

import unittest
import yaml
from pipeline_definition.pipeline_translator import PipelineTranslator

import examples.unix_commands


_yml = """
inputs:
  tar_input:
    tar:
      path: ../test-data/*.tar
      label: 'this is the input file'
"""

_expected_resolved = yaml.load("""
inputs:
  tar_input:
  - class: File
    path: ../test-data/hello.tar
""")

_expected_unresolved = yaml.load("""
inputs:
  tar_input:
    class: File
    path: ../test-data/*.tar
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
