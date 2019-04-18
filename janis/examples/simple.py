import unittest

from janis.hints import CaptureType

from janis.unix.tools.debugecho import DebugEchoInputs
from janis.unix.tools.tar import Tar
from janis.unix.tools.untar import Untar
from janis.unix.tools.compile import Compile
from janis.unix.data_types.tar_file import TarFile

from janis import Workflow, Input, Output, Step, String, File, Logger

# file.tar -> untar -> compile -> tar -> out.tar
#                  \_____________↗


class TestSimple(unittest.TestCase):

    def test_workflow(self):
        # Logger.mute()
        w = Workflow("simple")

        inp1 = Input("tarFile", TarFile())

        step1 = Step("untar", Untar())
        step2 = Step("compile", Compile())
        step3 = Step("tar", Tar())

        outp = Output("out", File())

        w.add_edge(inp1, step1.tarFile)
        w.add_edge(step1.files, step2.file)

        w.add_edge(inp1, step1.tarFile)
        w.add_edge(step1.files, step2.file)     # Auto scatter
        w.add_edge(step1.files, step3.files)
        w.add_edge(step2.compiled, step3.files)
        w.add_edge(step3.tarred, outp)
        # Logger.unmute()

        w.dump_translation("wdl", to_disk=False, with_resource_overrides=True)
        print(w.generate_resources_file("wdl", {CaptureType.KEY: CaptureType.CHROMOSOME}))


        return w

    def test_pipe(self):

        # Logger.mute()
        w = Workflow("simple")

        inp1 = Input("tarFile", TarFile(), "/Users/franklinmichael/Desktop/workflows-for-testing/03-simple/inputs/hello.tar")
        step1 = Step("untar", Untar())
        step2 = Step("compile", Compile())
        step3 = Step("tar", Tar())
        # debug = Step("debug", DebugEchoInputs())

        w.add_pipe(inp1, step1, step2, step3.files, Output("out"))
        w.add_edge(step1, step3.files2)
        # w.add_edge(step3.tarred, debug.fileInput)
        # w.add_edge(debug, Output("debugout"))

        # w.draw_graph()
        w.dump_translation(translation="cwl", to_disk=True, with_resource_overrides=False)
        # w.dump_wdl(to_disk=True, with_docker=True)

        # Logger.unmute()

        return w



