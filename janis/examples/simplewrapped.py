"""
simple.py

    Date: 2019-05-03
    Author: Michael Franklin


    This class represents a simple unix workflow, and was prepared as part of the
    following tutorial: https://janis.readthedocs.io/en/latest/tutorials/simple.html

    This version of the simple workflow IS wrapped in a subclass of Workflow to
    demonstrate how a subworkflow can be subclassed and easily be used again.

    Notes:
        - The toolshed only accepts workflows that are a subclass of `Workflow`,
        - It's easier to associate metadata with a subclass,
        - ** The runner (currently titled Shepherd) has a much easier time finding
                a subclass over a variable of type.


    Workflow description:
        1. A file is untarred
        2. The result of the untar is compiled
        3. The (untarred + compiled) result is tarred into a new archive

        The following ASCII diagram is a basic visual representation of the DAG:

        file.tar → untar → compile → tar -> out.tar
                          ↘_______↗
"""

# The classes we require to build a basic workflow
from janis import Workflow, Input, Output, Step

# Data types - These help us logically connect workflows
from janis.unix.data_types.tarfile import TarFile

# Tools - The command line tools we're going to call
from janis.unix.tools.compile import Compile
from janis.unix.tools.tar import Tar
from janis.unix.tools.untar import Untar


class SimpleWorkflow(Workflow):
    def __init__(self):
        super().__init__("simpleWorkflow")

        inp = Input("tarFile", TarFile(), value="/path/to/hello.tar")

        untar = Step("untar", Untar())
        compile = Step("compile", Compile())
        retar = Step("tar", Tar())

        outp = Output("out")

        self.add_edge(inp, untar.tarFile)
        self.add_edge(untar.out, compile.file)  # Auto scatter
        self.add_edge(untar.out, retar.files)
        self.add_edge(compile.out, retar.files2)
        self.add_edge(retar.out, outp)


if __name__ == "__main__":
    SimpleWorkflow().translate("wdl", to_disk=True)
