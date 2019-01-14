import unittest

from wdlgen.common import Input, Output
from wdlgen.task import Runtime, Command, Task
from wdlgen.types import WdlType, String
from wdlgen.workflow import Workflow
from wdlgen.workflowcall import WorkflowCall, WorkflowScatter


class TestTaskGeneration(unittest.TestCase):

    def test_simple_task(self):

        # Task based on: https://github.com/openwdl/wdl/blob/master/versions/draft-2/SPEC.md#introduction
        t = Task("hello", runtime=Runtime(docker="broadinstitute/my_image"))
        t.inputs.extend([
            Input(WdlType.parse_type("String"), "pattern"),
            Input(WdlType.parse_type("File"), "in")
        ])

        t.outputs.append(
            Output(WdlType.parse_type("Array[String]"), "matches", "read_lines(stdout())")
        )

        t.command = Command("egrep", inputs=[Command.CommandInput.from_input(i) for i in t.inputs])

        print(t.get_string())
        return t

    @staticmethod
    def test_hello_tasks():
        # https://github.com/openwdl/wdl/#specifying-inputs-and-using-declarations

        t1 = Task("hello",
                  inputs=[Input(String, "name"), Input(String, "salutation")],
                  outputs=[Output(String, "response", "read_string(stdout())")]
                  )
        t1.command = Command("echo", inputs=[Command.CommandInput.from_input(t1.inputs[i], position=i) for i in range(2)])

        print(t1.get_string())

        return t1


class TestCommandGeneration(unittest.TestCase):

    def test_simple_command(self):
        command = Command("egrep")
        command.inputs.append(Command.CommandInput("pattern"))
        command.inputs.append(Command.CommandInput("in"))

        print(command.get_string(2))


class TestWorkflowGeneration(unittest.TestCase):
    def test_hello_workflow(self):
        # https://github.com/openwdl/wdl/#specifying-inputs-and-using-declarations
        w = Workflow("test")

        w.calls.append(WorkflowCall(
            TestTaskGeneration.test_hello_tasks()
        ))

        w.calls.append(WorkflowCall(
            TestTaskGeneration.test_hello_tasks(),
            alias="hello2",
            inputs_map={
                "salutation": '"Greetings"',
                "name": '"Michael"'
            }
        ))

        print(w.get_string())


class TestWorkflowScatter(unittest.TestCase):
    def test_call_scatter(self):
        sc = WorkflowScatter("i", "integers", [
            WorkflowCall(Task("task1"), inputs_map={"num": "i"}),
            WorkflowCall(Task("task2"), inputs_map={"num": "task1.output"})
        ])

        print(sc.get_string(indent=0))
        print(sc.get_string(indent=1))
        print(sc.get_string(indent=2))
