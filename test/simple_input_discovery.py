
#
# Test discovery of inputs in the file system and generation of the input definition required by CWL.
#

import unittest
import yaml

from wehi.spec import Wehi

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
        translator = Wehi("simple_input_discovery")
        translator.parse_string(_yml)

        self.assertTrue(True)

    def test_resolved(self):
        self.translate(True, _expected_resolved)

    def test_unresolved(self):
        self.translate(False, _expected_unresolved)


if __name__ == '__main__':
    unittest.main()
