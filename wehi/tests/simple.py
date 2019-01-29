import unittest

from wehi.spec import Logger, Wehi


_simple = """
inputs:
    tarfile:
        type: TarFile
        path: /Users/franklinmichael/source/simple-workflow/hello.tar

steps:
    untar:
        tool: untar
        tarFile: tarfile    

    compile:
        tool: javacompiler
        file: untar     # will implicitly scatter

    tar:
        tool: tar
        files: compile  # will implicitly merge
        files2: untar 

outputs:
    untarred: untar
    compiled: compile
    tarred: tar
"""


class SimpleTest(unittest.TestCase):

    def test_simple(self, mute=False):
        import Pipeline.unix

        if mute:
            Logger.mute()

        w = Wehi("Simple")
        w.parse_string(_simple)
        w.workflow.dump_cwl(w.workflow)
        # w.workflow.dump_wdl(w.workflow, to_disk=True)
        Logger.unmute()
        self.assertTrue(True)
