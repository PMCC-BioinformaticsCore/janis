#
# This creates CWL for the first workflow example in in the CWL "documentation"
#

import yaml, unittest
from wehi.spec import Wehi

_yml = """
inputs:
  mytar:
    type: tar_file
    path: 'test-data/hello.tar'

steps:
  untar_step:
    tool: untar
    input: mytar
    
  compile_step:
    tool: java-compile
    input: 'untar_step/out'
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
    translator = Wehi("1st-workflow")
    translator.parse_string(_yml)
    print('-'*80)
    print(translator.workflow.cwl()[0])
    print('-'*80)
    self.assertTrue(True)


if __name__ == '__main__':
    FirstWorkflow().test_graph()

