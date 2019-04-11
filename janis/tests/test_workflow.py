from unittest import TestCase

from janis import File, Array, Logger, CommandTool, ToolInput, String, Input, Output, Step, Workflow
from janis.graph.stepinput import StepInput, first_value, Edge

from janis.unix.data_types.tar_file import TarFile
from janis.unix.tools.echo import Echo
from janis.unix.tools.tar import Tar
from janis.unix.tools.untar import Untar
from janis.unix.tools.cat import Cat


class SingleTestTool(CommandTool):

    @staticmethod
    def tool():
        return "TestStepTool"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [
            ToolInput("inputs", String())
        ]

    def friendly_name(self):
        return None

    def outputs(self):
        return []

    @staticmethod
    def docker():
        return None


class ArrayTestTool(CommandTool):

    @staticmethod
    def tool():
        return "ArrayStepTool"

    def friendly_name(self):
        return None

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [
            ToolInput("inputs", Array(String()))
        ]

    def outputs(self):
        return []

    @staticmethod
    def docker():
        return None


class TestWorkflow(TestCase):

    def setUp(self):
        Logger.mute()

    def tearDown(self):
        Logger.unmute()

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
        w._add_item(inp)
        self.assertEqual(len(w._inputs), 1)
        self.assertEqual(w._inputs[0].input, inp)
        self.assertIsNotNone(w._nodes[inp.id()])

    def test_add_step(self):
        w = Workflow("test_add_input")
        step = Step("catStep", Cat())
        w._add_item(step)
        self.assertEqual(len(w._steps), 1)
        self.assertEqual(w._steps[0].step, step)
        self.assertIsNotNone(w._nodes[step.id()])

    def test_add_output(self):
        w = Workflow("test_add_input")
        outp = Output("outputStep", String())
        w._add_item(outp)
        self.assertEqual(len(w._outputs), 1)
        self.assertEqual(w._outputs[0].output, outp)
        self.assertIsNotNone(w._nodes[outp.id()])

    def test_add_node(self):
        w = Workflow("test_add_node")
        inp = Input("inp", String())
        stp = Step("stp", Cat())
        w.add_items([inp, stp])
        self.assertEqual(len(w._inputs), 1)
        self.assertEqual(len(w._steps), 1)
        self.assertEqual(len(w._outputs), 0)
        self.assertEqual(w._nodes[inp.id()].id(), inp.id())
        self.assertEqual(w._nodes[stp.id()].id(), stp.id())

    def test_add_qualified_edge(self):
        w = Workflow("test_add_edge")
        inp = Input("inp", String())
        stp = Step("stp", Echo())  # Only has one input, with no output
        e = w.add_edge(inp, stp.input)

        self.assertEqual(e.start.id(), inp.id())
        self.assertEqual(e.finish.id(), stp.id())
        self.assertIsNone(e.stag)
        self.assertEqual(e.ftag, next(iter(stp.requires())))

    def test_add_edge(self):
        w = Workflow("test_add_edge")
        inp = Input("inp", String())
        stp = Step("stp", Echo())       # Only has one input, with no output
        w.add_items([inp, stp])
        e = w.add_edge(inp, stp)

        self.assertEqual(e.start.id(), inp.id())
        self.assertEqual(e.finish.id(), stp.id())
        self.assertIsNone(e.stag)
        self.assertEqual(e.ftag, next(iter(stp.requires())))

    def test_anonymous_add_edge(self):
        w = Workflow("test_add_edge")
        inp = Input("inp", String())
        stp = Step("stp", Echo())       # Only has one input, with no output
        # w.add_items([inp, stp])
        e = w.add_edge(inp, stp)

        self.assertEqual(e.start.id(), inp.id())
        self.assertEqual(e.finish.id(), stp.id())
        self.assertIsNone(e.stag)
        self.assertEqual(e.ftag, next(iter(stp.requires())))

    def test_anonymous_add_qualified_edge(self):
        w = Workflow("test_add_edge")
        inp = Input("inp", String())
        stp = Step("stp", Echo())       # Only has one input, with no output
        e = w.add_edge(inp, stp.input)

        self.assertEqual(e.start.id(), inp.id())
        self.assertEqual(e.finish.id(), stp.id())
        self.assertIsNone(e.stag)
        self.assertEqual(e.ftag, next(iter(stp.requires())))

    def test_pipe(self):
        w = Workflow("test_add_edge")
        inp = Input("tarred", TarFile())
        stp = Step("stp", Untar())  # Only has one input, with no output
        out = Output("outp", Array(File()))

        w.add_pipe(inp, stp, out)

        # the nodes are usually internal
        inp_node = w._nodes[inp.id()]
        stp_node = w._nodes[stp.id()]
        out_node = w._nodes[out.id()]

        self.assertEqual(len(inp_node.connection_map), 0)
        self.assertEqual(len(stp_node.connection_map), 1)
        self.assertEqual(len(out_node.connection_map), 1)

        s1: StepInput = first_value(stp_node.connection_map)
        s2: StepInput = first_value(out_node.connection_map)

        e1 = first_value(s1.source_map)
        e2 = first_value(s2.source_map)

        self.assertEqual(e1.start.id(),  inp.id())
        self.assertEqual(e1.finish.id(), stp.id())
        self.assertEqual(e2.start.id(),  stp.id())
        self.assertEqual(e2.finish.id(), out.id())

    def test_qualified_pipe(self):
        w = Workflow("test_add_edge")
        inp = Input("tarred", TarFile())
        stp = Step("stp", Untar())  # Only has one input, with no output
        out = Output("outp", Array(File()))

        w.add_pipe(inp, stp.files, out)

        # the nodes are usually internal
        inp_node = w._nodes[inp.id()]
        stp_node = w._nodes[stp.id()]
        out_node = w._nodes[out.id()]

        self.assertEqual(len(inp_node.connection_map), 0)
        self.assertEqual(len(stp_node.connection_map), 1)
        self.assertEqual(len(out_node.connection_map), 1)

        s1: StepInput = first_value(stp_node.connection_map)
        s2: StepInput = first_value(out_node.connection_map)

        e1: Edge = first_value(s1.source_map)
        e2: Edge = first_value(s2.source_map)

        self.assertEqual(e1.start.id(), inp.id())
        self.assertEqual(e1.finish.id(), stp.id())
        self.assertEqual(e2.start.id(), stp.id())
        self.assertEqual(e2.finish.id(), out.id())

    def test_subworkflow(self):

        w = Workflow("test_subworkflow")

        sub_w = Workflow("subworkflow")
        sub_inp = Input("sub_inp", TarFile())
        sub_stp = Step("sub_stp", Untar())
        sub_out = Output("sub_out", Array(File()))
        sub_w.add_pipe(sub_inp, sub_stp, sub_out)

        inp = Input("inp", TarFile())
        stp = Step("stp_workflow", sub_w)
        out = Output("out", Array(File()))
        w.add_items([inp, stp, out])
        w.add_pipe(inp, stp, out)

        # w.dump_cwl(to_disk=True)

        self.assertTrue(True)


    def test_add_scatter(self):
        w = Workflow("scatterededge")

        inp1 = Input("inp1", Array(String()))
        step = Step("stp", SingleTestTool())

        e = w.add_edge(inp1, step)
        self.assertTrue(e.scatter)

    def test_add_scatter_nested_arrays(self):
        w = Workflow("scatterededge")

        inp1 = Input("inp1", Array(Array(String())))
        step = Step("stp", ArrayTestTool())

        e = w.add_edge(inp1, step)
        self.assertTrue(e.scatter)

    def test_add_non_scatter(self):
        w = Workflow("scatterededge")

        inp1 = Input("inp1", String())
        step = Step("stp", SingleTestTool())

        e = w.add_edge(inp1, step)
        self.assertFalse(e.scatter)

    def test_add_non_scatter2(self):
        w = Workflow("scatterededge")

        inp1 = Input("inp1", Array(String()))
        step = Step("stp", ArrayTestTool())

        e = w.add_edge(inp1, step)
        self.assertFalse(e.scatter)
