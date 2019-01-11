import unittest

import wdlgen.task as task
from wdlgen.types import WdlType


class TestTaskGeneration(unittest.TestCase):

    def test_simple_task(self):
        t = task.Task("hello", runtime=task.Runtime(docker="broadinstitute/my_image"))
        t.inputs.extend([
            task.Input(WdlType.parse_type("String"), "pattern"),
            task.Input(WdlType.parse_type("File"), "in")
        ])

        t.outputs.append(
            task.Output(WdlType.parse_type("Array[String]"), "matches", "read_lines(stdout())")
        )

        t.command = task.Command("egrep", inputs=[task.Command.CommandInput.from_input(i) for i in t.inputs])

        print(t.wdl())


class TestCommandGeneration(unittest.TestCase):

    def test_simple_command(self):
        command = task.Command("egrep")
        command.inputs.append(task.Command.CommandInput("pattern"))
        command.inputs.append(task.Command.CommandInput("in"))

        print(command.wdl(2))

