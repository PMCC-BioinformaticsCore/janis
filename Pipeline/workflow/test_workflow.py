from unittest import TestCase

from Pipeline import File, Array
from Pipeline.types.common_data_types import String
from Pipeline.unix.data_types.tar_file import TarFile
from Pipeline.unix.steps.echo import Echo
from Pipeline.unix.steps.untar import Untar
from Pipeline.workflow.input import Input
from Pipeline.workflow.output import Output
from Pipeline.workflow.step import Step
from Pipeline.workflow.workflow import Workflow
from Pipeline.unix.steps.cat import Cat


class TestWorkflow(TestCase):

    def test_name(self):
        wn = "test_name"
        w = Workflow(wn)
        self.assertEqual(w.identifier, wn)

    def test_rename(self):
        wn1 = "test_rename"
        wn2 = "test_rename2"
        w = Workflow(wn1)
        w.identifier = wn2
        self.assertEqual(w.identifier, wn2)

    def test_add_input(self):
        w = Workflow("test_add_input")
        inp = Input("inputLabel", String())
        w.add_input(inp)
        self.assertEqual(len(w._inputs), 1)
        self.assertEqual(w._inputs[0].input, inp)
        self.assertIsNotNone(w._labels[inp.id()])

    def test_add_step(self):
        w = Workflow("test_add_input")
        step = Step("catStep", Cat())
        w.add_step(step)
        self.assertEqual(len(w._steps), 1)
        self.assertEqual(w._steps[0].step, step)
        self.assertIsNotNone(w._labels[step.id()])

    def test_add_output(self):
        w = Workflow("test_add_input")
        outp = Output("outputStep", String())
        w.add_output(outp)
        self.assertEqual(len(w._outputs), 1)
        self.assertEqual(w._outputs[0].output, outp)
        self.assertIsNotNone(w._labels[outp.id()])

    def test_add_node(self):
        w = Workflow("test_add_node")
        inp = Input("inp", String())
        stp = Step("stp", Cat())
        w.add_nodes([inp, stp])
        self.assertEqual(len(w._inputs), 1)
        self.assertEqual(len(w._steps), 1)
        self.assertEqual(len(w._outputs), 0)
        self.assertEqual(w._labels[inp.id()].id(), inp.id())
        self.assertEqual(w._labels[stp.id()].id(), stp.id())

    def test_add_qualified_edge(self):
        w = Workflow("test_add_edge")
        inp = Input("inp", String())
        stp = Step("stp", Echo())  # Only has one input, with no output
        w.add_nodes([inp, stp])
        e = w.add_edge(inp, stp.inp)

        self.assertEqual(e.start.id(), inp.id())
        self.assertEqual(e.finish.id(), stp.id())
        self.assertIsNone(e.stag)
        self.assertEqual(e.ftag, next(iter(stp.requires())))

    def test_add_edge(self):
        w = Workflow("test_add_edge")
        inp = Input("inp", String())
        stp = Step("stp", Echo())       # Only has one input, with no output
        w.add_nodes([inp, stp])
        e = w.add_edge(inp, stp)

        self.assertEqual(e.start.id(), inp.id())
        self.assertEqual(e.finish.id(), stp.id())
        self.assertIsNone(e.stag)
        self.assertEqual(e.ftag, next(iter(stp.requires())))

    def test_pipe(self):
        w = Workflow("test_add_edge")
        inp = Input("tarred", TarFile())
        stp = Step("stp", Untar())  # Only has one input, with no output
        outp = Output("outp", Array(File()))

        w.add_nodes([inp, stp, outp])
        w.add_pipe(inp, stp, outp)

        # the nodes are usually internal
        inp_node = w._labels[inp.id()]
        stp_node = w._labels[stp.id()]
        out_node = w._labels[outp.id()]



    def test_qualified_pipe(self):
        w = Workflow("test_add_edge")
        inp = Input("tarred", TarFile())
        stp = Step("stp", Untar())  # Only has one input, with no output
        outp = Output("outp", Array(File()))

        w.add_nodes([inp, stp, outp])
        w.add_pipe(inp, stp.files, outp)





