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
from janis_core import Workflow, Array

# Data types - These help us logically connect workflows
from janis.unix.data_types.tarfile import TarFile

# Tools - The command line tools we're going to call
from janis.unix.tools.compile import Compile
from janis.unix.tools.tar import Tar
from janis.unix.tools.untar import Untar


w = Workflow("simple")

w.input(
    "tarFile",
    TarFile,
    default="/Users/franklinmichael/Desktop/workflows-for-testing/03-simple/inputs/hello.tar",
)

w.step("untar", Untar, tarFile=w.tarFile)
w.step("compile", Compile, scatter="file", file=w.untar.out)
w.step("tar", Tar, files=w.untar.out, files2=w.compile.out)

w.output("outp", source=w.tar.out)

if __name__ == "__main__":
    w.translate("wdl", to_disk=True, validate=True)
