from janis import Workflow, Input, File, Output, Step
from janis.unix.data_types.tar_file import TarFile
from janis.unix.tools.compile import Compile
from janis.unix.tools.untar import Untar

w = Workflow("user_guide")

tarball = Input("tarball", TarFile())
compiled_class = Output("compiled_class")

untar = Step("untar", Untar())
compile_step = Step("compile", Compile())

w.add_edge(tarball, untar.tarFile)
w.add_edge(untar, compile_step)
w.add_edge(compile_step, compiled_class)

w.dump_translation("cwl", to_disk=False)
w.dump_wdl(to_disk=True)