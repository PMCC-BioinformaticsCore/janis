"""
simple.py

    Date: 2019-05-03
    Author: Michael Franklin


    This file contains a simple unix workflow, and was prepared as part of the
    following tutorial: https://janis.readthedocs.io/en/latest/tutorials/simple.html

    This version of the simple workflow is NOT wrapped in a subclass of Workflow.

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


w = Workflow("simple")

inp = Input("tarFile", TarFile(), value="/path/to/hello.tar")

untar = Step("untar", Untar())
compil = Step("compile", Compile())
retar = Step("tar", Tar())

outp = Output("out")

w.add_edge(inp, untar.tarFile)
w.add_edge(untar.files, compil.file)  # Auto scatter
w.add_edge(untar.files, retar.files)
w.add_edge(compil.compiled, retar.files2)
w.add_edge(retar.tarred, outp)


if __name__ == "__main__":
    w.translate("wdl", to_disk=True)
