import unittest

from wdlgen import Input, Output, Workflow, Task, WdlType, String, WorkflowCall, WorkflowScatter


class TestTaskGeneration(unittest.TestCase):

    def test_simple_task(self):
        # Task based on: https://github.com/openwdl/wdl/blob/master/versions/draft-2/SPEC.md#introduction
        t = Task("hello", runtime=Task.Runtime(docker="broadinstitute/my_image"))
        t.inputs.extend([
            Input(WdlType.parse_type("String"), "pattern"),
            Input(WdlType.parse_type("File"), "in")
        ])

        t.outputs.append(
            Output(WdlType.parse_type("Array[String]"), "matches", "read_lines(stdout())")
        )

        t.command = Task.Command("egrep", inputs=[Task.Command.CommandInput.from_input(i) for i in t.inputs])

        print(t.get_string())
        return t

    def test_readme(self):
        w = Workflow("workflow_name")

        w.imports.append(Workflow.WorkflowImport("tool_file", ""))
        w.inputs.append(Input(WdlType.parse_type("String"), "inputGreeting"))

        inputs_map = {"taskGreeting": "inputGreeting"}
        w.calls.append(WorkflowCall("Q.namspaced_task_identifier", "task_alias", inputs_map))
        w.outputs.append(Output(WdlType.parse_type("File"), "standardOut", "task_alias.standardOut"))

        print(w.get_string())

    def test_readme_task(self):
        t = Task("task_name")
        t.inputs.append(Input(WdlType.parse_type("String"), "taskGreeting"))
        # command in next section
        t.outputs.append(Output(WdlType.parse_type("File"), "standardOut", "stdout()"))

        command = Task.Command("echo")
        command.inputs.append(
            Task.Command.CommandInput("taskGreeting", optional=False, position=None, prefix="-a",
                                             separate_value_from_prefix=True, default=None))
        command.inputs.append(
            Task.Command.CommandInput("otherInput", optional=True, position=2, prefix="optional-param=",
                                             separate_value_from_prefix=False, default=None))
        command = Task.Command("echo")
        command.inputs.append(
            Task.Command.CommandInput("taskGreeting", optional=False, position=None, prefix="-a",
                                      separate_value_from_prefix=True, default=None))
        command.inputs.append(
            Task.Command.CommandInput("otherInput", optional=True, position=2, prefix="optional-param=",
                                      separate_value_from_prefix=False, default=None))

        # t is the task
        t.command = command
        print(t.get_string())


    @staticmethod
    def test_hello_tasks():
        # https://github.com/openwdl/wdl/#specifying-inputs-and-using-declarations

        t1 = Task("hello",
                  inputs=[Input(String, "name"), Input(String, "salutation")],
                  outputs=[Output(String, "response", "read_string(stdout())")]
                  )
        t1.command = Task.Command("echo",
                                  inputs=[Task.Command.CommandInput.from_input(t1.inputs[i], position=i) for i in
                                          range(2)])

        print(t1.get_string())

        return t1


class TestCommandGeneration(unittest.TestCase):

    def test_simple_command(self):
        command = Task.Command("egrep")
        command.inputs.append(Task.Command.CommandInput("pattern"))
        command.inputs.append(Task.Command.CommandInput("in"))

        print(command.get_string(2))

    def test_readme_example(self):
        command = Task.Command("echo")
        command.inputs.append(
            Task.Command.CommandInput("taskGreeting", optional=False, position=None, prefix="-a",
                                      separate_value_from_prefix=True, default=None))
        command.inputs.append(
            Task.Command.CommandInput("otherInput", optional=True, position=2, prefix="optional-param=",
                                      separate_value_from_prefix=False, default=None))

        # t is the task
        print(command.get_string())

    def test_commandinput_space(self):
        t = Task.Command.CommandInput("taskGreeting", optional=False, position=None, prefix="-a",
                                      separate_value_from_prefix=True, default=None)
        self.assertEqual("-a ${taskGreeting}", t.get_string())

    def test_commandinput_nospace(self):
        t = Task.Command.CommandInput("taskGreeting", optional=False, position=None, prefix="val=",
                                      separate_value_from_prefix=False, default=None)
        self.assertEqual("val=${taskGreeting}", t.get_string())

    def test_commandarg_space(self):
        t = Task.Command.CommandInput("argVal", position=None, prefix="-p", separate_value_from_prefix=True)
        self.assertEqual("-p ${argVal}", t.get_string())

    def test_commandarg_nospace(self):
        t = Task.Command.CommandArgument(prefix="arg=", value="argVal", position=None, separate_value_from_prefix=False)
        self.assertEqual("arg=argVal", t.get_string())


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
