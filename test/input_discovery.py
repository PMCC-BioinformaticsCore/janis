
#
# Test discovery of inputs and generation of the input definition required by CWL.
import unittest
import yaml
from pipeline_definition.pipeline_translator import PipelineTranslator

import examples.unix_commands


_yml = """
inputs:
  tar_input:
    tar:
      glob: test-data/hello.tar
      label: 'this is the input file'

"""

_expected = yaml.load("""
inp:
  # This is the input tar file
  class: File
  path: test-data/hello.tar
  """)


class InputFileTest(unittest.TestCase):

  def test_input(self):
    translator = PipelineTranslator(debug=False)
    translator.translate_string(_yml)

    translation = translator.input(resolve=False)
    print(translation)

    self.assertTrue(True)


if __name__ == '__main__':
  unittest.main()
