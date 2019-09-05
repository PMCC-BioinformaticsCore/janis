from janis import Workflow
from janis.unix.data_types.tarfile import TarFile
from janis.unix.tools.compile import Compile
from janis.unix.tools.untar import Untar

w = Workflow("user_guide")

tarball = Input("tarball", TarFile(), value="path/to/tar")
compiled_class = Output("compiledOutput")

untar = Step("untar", Untar())
compile_step = Step("compile", Compile())

w.add_edges(
    [
        (tarball, untar.tarFile),
        (untar.out, compile_step.file),
        (compile_step.out, compiled_class),
    ]
)

w.translate("cwl", to_disk=False)  # or "wdl"
