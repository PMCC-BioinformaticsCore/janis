from unittest import TestCase

from Pipeline.types.common_data_types import String
from Pipeline.workflow.input import Input
from Pipeline.workflow.output import Output
from Pipeline.workflow.step import Step
from Pipeline.workflow.workflow import Workflow
from examples.unix_commands.steps.cat import Cat


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

