import unittest

from janis.unix.tools.debugecho import DebugEchoInputs
from janis.unix.tools.tar import Tar
from janis.unix.tools.untar import Untar
from janis.unix.tools.compile import Compile
from janis.unix.data_types.tar_file import TarFile

from janis import Workflow, Input, Output, Step, String, File, Logger

# file.tar -> untar -> compile -> tar -> out.tar
#                  \_____________↗


class TestSimple(unittest.TestCase):

    # def test_workflow(self):
    #     Logger.mute()
    #     w = Workflow("simple")
    #
    #     inp1 = Input("tarFile", TarFile())
    #     # inp2 = Input("tarName", String())
    #
    #     step1 = Step("untar", Untar())
    #     step2 = Step("compile", Compile())
    #     step3 = Step("tar", Tar())
    #
    #     outp = Output("output", File())
    #
    #     w.add_edge(inp1, step1.tarFile)
    #     w.add_edge(step1.files, step2.file)
    #
    #     w.add_edge(inp1, step1.tarFile)
    #     # w.add_edge(inp2, step3.tarName)
    #     w.add_edge(step1.files, step2.file)     # Auto scatter
    #     w.add_edge(step1.files, step3.files)
    #     w.add_edge(step2.compiled, step3.files)
    #     w.add_edge(step3.tarred, outp)
    #     Logger.unmute()
    #
    #     w.dump_wdl(to_disk=False)
    #
    #     return w

    def test_pipe(self):

        # Logger.mute()
        w = Workflow("simple")

        inp1 = Input("tarFile", TarFile(), "/Users/franklinmichael/source/simple-workflow/hello.tar")
        step1 = Step("untar", Untar())
        step2 = Step("compile", Compile())
        step3 = Step("tar", Tar())
        debug = Step("debug", DebugEchoInputs())

        w.add_pipe(inp1, step1, step2, step3.files, Output("output"))
        w.add_edge(step1, step3.files2)
        w.add_edge(step3, debug)
        w.add_edge(debug, Output("debugout"))

        # w.draw_graph()
        w.dump_translation(translation="cwl", to_disk=False, with_docker=False)
        # w.dump_wdl(to_disk=True, with_docker=True)

        # Logger.unmute()

        return w



