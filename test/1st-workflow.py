#
# This creates CWL for the first workflow example in in the CWL "documentation"
#

import unittest
import yaml
from pipeline_definition.pipeline_translator import PipelineTranslator
import examples.unix_commands


_yml = """
inputs:
  mytar:
    tar_file:
      path: 'test-data/hello.tar'

steps:
  - untar_step:
      untar:
  - compile_step:
      compile
"""

_expected_cwl = yaml.load("""
cwlVersion: v1.0
class: Workflow
inputs:
  inp: File
  ex: string

outputs:
  classout:
    type: File
    outputSource: compile/classfile

steps:
  untar:
    run: tar-param.cwl
    in:
      tarfile: inp
      extractfile: ex
    out: [example_out]

  compile:
    run: arguments.cwl
    in:
      src: untar/example_out
    out: [classfile]
""")

_expect_inp = yaml.load("""
inp:
  class: File
  path: test-data/hello.tar
ex: Hello.java
""")


class FirstWorkflow(unittest.TestCase):

  def test_graph(self):
    translator = PipelineTranslator(debug=True)
    translation = translator.translate_string(_yml)

    print('-'*80)
    print(translation)
    print('-'*80)

    self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()

