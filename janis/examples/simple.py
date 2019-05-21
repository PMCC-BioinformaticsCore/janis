"""
simple.py

    Date: 2019-05-03
    Author: Michael Franklin


    This class represents a simple unix workflow, and was prepared as part of the
    following tutorial: https://janis.readthedocs.io/en/latest/tutorials/simple.html

    Workflow description:
        1. A file is untarred
        2. The result of the untar is compiled
        3. The (untarred + compiled) result is tarred into a new archive

        The following ASCII diagram is a basic visual representation of the DAG:

        file.tar → untar → compile → tar -> out.tar
                          ↘_______↗
"""

# The classes we require to build a basic workflow
from janis import Workflow, Input, Output, Step, File

# Data types - These help us logically connect workflows
from janis.unix.data_types.tarfile import TarFile

# Tools - The command line tools we're going to call
from janis.unix.tools.compile import Compile
from janis.unix.tools.tar import Tar
from janis.unix.tools.untar import Untar


class SimpleWorkflow(Workflow):

    def __init__(self):
        super().__init__("simpleWorkflow")

        inp = Input("tarFile", TarFile())

        step1 = Step("untar", Untar())
        step2 = Step("compile", Compile())
        step3 = Step("tar", Tar())

        outp = Output("out")

        self.add_edge(inp, step1.tarFile)
        self.add_edge(step1.files, step2.file)  # Auto scatter
        self.add_edge(step1.files, step3.files)
        self.add_edge(step2.compiled, step3.files2)
        self.add_edge(step3.tarred, outp)


if __name__ == "__main__":
    SimpleWorkflow().translate("wdl", to_disk=True, merge_resources=True)



