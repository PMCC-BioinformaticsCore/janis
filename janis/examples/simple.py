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
from janis import Workflow, Input, Output, Step, Array

# Data types - These help us logically connect workflows
from janis.unix.data_types.tarfile import TarFile

# Tools - The command line tools we're going to call
from janis.unix.tools.compile import Compile
from janis.unix.tools.tar import Tar
from janis.unix.tools.untar import Untar


w = Workflow("simple")

inp = Input(
    "tarFile",
    TarFile(),
    default="/Users/franklinmichael/Desktop/workflows-for-testing/03-simple/inputs/hello.tar",
)

untar = Step("untar", Untar())
compil = Step("compile", Compile())
tar = Step("tar", Tar())

outp = Output("out")

w.add_edge(inp, untar.tarFile)
w.add_edge(untar.out, compil.file)  # Auto scatter
w.add_edge(untar.out, tar.files)
w.add_edge(compil.out, tar.files2)
w.add_edge(tar.out, outp)


if __name__ == "__main__":
    w.translate("wdl", to_disk=True, should_validate=True)
