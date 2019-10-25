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
from janis_core import Workflow

# Data types - These help us logically connect workflows
from janis.data_types import TarFile

# Tools - The command line tools we're going to call
from janis.tools import Compile, Untar


class SimpleWorkflow(Workflow):
    def id(self) -> str:
        return "simple"

    def friendly_name(self):
        return "Simple workflow: Untar, Compile, Tar"

    def constructor(self):
        self.input("tarfile", TarFile)
        self.step("untar", Untar(tarfile=self.tarfile))
        self.step("compile", Compile(file=self.untar.out), scatter="file")
        self.output("out", source=self.compile)


if __name__ == "__main__":
    SimpleWorkflow().translate("wdl", to_disk=True)
