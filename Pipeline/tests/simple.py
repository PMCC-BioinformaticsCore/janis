import unittest

from Pipeline import Workflow, Input, Output, Step, String, File

from examples.unix_commands import Tar, Compile, Untar
from examples.unix_commands import TarFile

# Write simple workflow here


# file.tar -> untar -> compile -> tar -> out.tar
#                  \_____________↗


class TestSimple(unittest.TestCase):

    def test_workflow(self):
        w = Workflow("simple")

        inp1 = Input("tarFile", TarFile())
        inp2 = Input("tarName", String())

        step1 = Step("untar", Untar())
        step2 = Step("compile", Compile())
        step3 = Step("tar", Tar())

        outp = Output("output", File())

        w.add_input(inp1)
        w.add_input(inp2)
        w.add_step(step1)
        w.add_step(step2)
        w.add_step(step3)
        w.add_output(outp)

        w.add_edge(inp1, step1)
        w.add_edge(step1, step2)
        # w.add_pipe([inp1, step1, step2, step3])

        w.add_edge(inp1, step1.tarFile)
        w.add_edge(inp2, step3.tarName)
        w.add_edge(step1.outp, step2.file)
        w.add_edge(step1.outp, step3.input1)
        w.add_edge(step2.outp, step3.input2)
        w.add_edge(step3.outp, outp)

        print(w.cwl())

        # w.draw_graph()
        return w
