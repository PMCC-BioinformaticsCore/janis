"""
simple.py

    Date: 2019-10-03
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
from janis import WorkflowBuilder
from janis.data_types import TarFile
from janis.tools import Compile, Untar

w = WorkflowBuilder("user_guide")

w.input("tarfile", TarFile, value="/path/to/file.tar")
w.step("untar", Untar(tarfile=w.tarfile))
w.step("compile", Compile(file=w.untar), scatter="file")
w.output("compiled", source=w.compile.out)

if __name__ == "__main__":
    w.translate("wdl", to_disk=True, validate=True)
