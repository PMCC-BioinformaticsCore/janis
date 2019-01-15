import unittest

from Pipeline.unix.steps.tar import Tar
from Pipeline.unix.steps.untar import Untar
from Pipeline.unix.steps.compile import Compile
from Pipeline.unix.data_types.tar_file import TarFile

from Pipeline import Workflow, Input, Output, Step, String, File, Logger


# file.tar -> untar -> compile -> tar -> out.tar
#                  \_____________↗


class TestSimple(unittest.TestCase):

    def test_workflow(self):
        Logger.mute()
        w = Workflow("simple")

        inp1 = Input("tarFile", TarFile())
        # inp2 = Input("tarName", String())

        step1 = Step("untar", Untar())
        step2 = Step("compile", Compile())
        step3 = Step("tar", Tar())

        outp = Output("output", File())

        w.add_edge(inp1, step1.tarFile)
        w.add_edge(step1.files, step2.file)

        w.add_edge(inp1, step1.tarFile)
        # w.add_edge(inp2, step3.tarName)
        w.add_edge(step1.files, step2.file)     # Auto scatter
        w.add_edge(step1.files, step3.files)
        w.add_edge(step2.compiled, step3.files)
        w.add_edge(step3.tarred, outp)
        Logger.unmute()

        w.dump_wdl(to_disk=False)

        return w

    def test_pipe(self):
        Logger.mute()
        w = Workflow("simple")

        inp1 = Input("tarFile", TarFile(), "/Users/franklinmichael/source/simple-workflow/hello.tar")
        step1 = Step("untar", Untar())
        step2 = Step("compile", Compile())
        step3 = Step("tar", Tar())
        outp = Output("output", File())

        w.add_pipe(inp1, step1, step2, step3.files, outp)
        w.add_edge(step1, step3.files)

        # w.draw_graph()
        w.dump_cwl(to_disk=True)
        w.dump_wdl(to_disk=True)

        Logger.unmute()

        return w



