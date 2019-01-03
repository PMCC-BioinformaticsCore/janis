from unittest import TestCase

from Pipeline import Input, String, CommandTool, ToolInput, Step, Array, Logger
from Pipeline.graph.stepinput import StepInput
from Pipeline.workflow.input import InputNode
from Pipeline.workflow.step import StepNode


class TestTool(CommandTool):

    @staticmethod
    def tool():
        return "TestStepTool"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [
            ToolInput("inputs", Array(String()))
        ]

    def outputs(self):
        return []


class TestStep(TestCase):

    def setUp(self):
        Logger.mute()

    def tearDown(self):
        Logger.unmute()

    def test_sources_single(self):
        start1 = InputNode(Input("test1", String()))
        step = StepNode(Step("step", TestTool()))

        con_map = StepInput(step, "inputs")
        con_map.add_source(start1, None)

        self.assertEqual(con_map.source(), start1.id())

    def test_sources_multiple(self):
        start1 = InputNode(Input("test1", String()))
        start2 = InputNode(Input("test2", String()))
        step = StepNode(Step("step", TestTool()))

        con_map = StepInput(step, "inputs")
        con_map.add_source(start1, None)
        con_map.add_source(start2, None)

        self.assertIsInstance(con_map.source(), list)
        self.assertEqual(len(con_map.source()), 2)

    def test_sources_none(self):
        step = StepNode(Step("step", TestTool()))
        con_map = StepInput(step, "inputs")

        self.assertIsNone(con_map.source())
