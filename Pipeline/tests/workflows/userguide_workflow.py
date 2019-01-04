from Pipeline import Workflow, Input, File, Output, Step
from Pipeline.unix.data_types.tar_file import TarFile
from Pipeline.unix.steps.compile import Compile
from Pipeline.unix.steps.tar import Tar
from Pipeline.unix.steps.untar import Untar

w = Workflow("user-guide")

tarball = Input("tarball", TarFile())
compiled_class = Output("compiled_class")

untar = Step("untar", Untar())
compile_step = Step("compile", Compile())

w.add_edge(tarball, untar)
w.add_edge(untar, compile_step)
w.add_edge(compile_step, compiled_class)

w.dump_cwl(to_disk=False)