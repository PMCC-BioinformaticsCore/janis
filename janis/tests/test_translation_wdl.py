import unittest
import wdlgen
from typing import List

from janis.translations import WdlTranslator

from janis.types import CpuSelector, MemorySelector, StringFormatter

import janis.translations.wdl as wdl

from janis import Workflow, Input, ToolOutput, ToolInput, String, CommandTool, Stdout, InputSelector, Array, File, Filename, \
    WildcardSelector


class TestTool(CommandTool):

    @staticmethod
    def tool(): return "TestTranslation-tool"

    @staticmethod
    def base_command(): return "echo"

    def inputs(self) -> List[ToolInput]: return [ToolInput("testtool", String())]

    def outputs(self) -> List[ToolOutput]: return [ToolOutput("std", Stdout())]

    def friendly_name(self) -> str: return "Tool for testing translation"

    @staticmethod
    def docker(): return "ubuntu:latest"


class TestToolWithSecondaryOutput(TestTool):
    def outputs(self):
        return [
            ToolOutput("out", TestTypeWithNonEscapedSecondary(), glob=InputSelector("testtool") + "/out")
        ]


class TestTypeWithSecondary(File):
    @staticmethod
    def secondary_files():
        return ["^.txt"]


class TestTypeWithNonEscapedSecondary(File):
    @staticmethod
    def secondary_files():
        return [".txt"]


class TestWdl(unittest.TestCase):

    def test_optional_array(self):
        t = Array(File(), optional=True)
        wdl = t.wdl()
        self.assertIsInstance(wdl, wdlgen.WdlType)
        self.assertTrue(wdl.optional)
        self.assertEqual("Array[File]?", wdl.get_string())


class TestWdlTranslatorOverrides(unittest.TestCase):

    def setUp(self):
        self.translator = WdlTranslator()

    def test_stringify_workflow(self):
        wdlobj = wdlgen.Workflow("wid", version="development")
        self.assertEqual(
            "version development\n\n\n\nworkflow wid {\n\n\n\n}",
            self.translator.stringify_translated_workflow(wdlobj)
        )

    def test_stringify_tool(self):
        wdlobj = wdlgen.Task("tid", version="development")
        self.assertEqual(
            "version development\n\ntask tid {\n\n\n\n\n}",
            self.translator.stringify_translated_tool(wdlobj)
        )

    def test_stringify_inputs(self):
        d = {"wid.inp1": 1}
        self.assertEqual(
            "{\n    \"wid.inp1\": 1\n}",
            self.translator.stringify_translated_inputs(d)
        )

    def test_workflow_filename(self):
        w = Workflow("wid")
        self.assertEqual("wid.wdl", self.translator.workflow_filename(w))

    def test_tools_filename(self):
        self.assertEqual("TestTranslation-tool.wdl", self.translator.tool_filename(TestTool().id()))

    def test_inputs_filename(self):
        w = Workflow("wid")
        self.assertEqual("wid-inp.json", self.translator.inputs_filename(w))

    def test_resources_filename(self):
        w = Workflow("wid")
        self.assertEqual("wid-resources.json", self.translator.resources_filename(w))


class TestWdlTranslatorBuilders(unittest.TestCase):

    def test_inputs_generator_secondary_files(self):
        w = Workflow("tst")
        w._add_input(Input("wsec", TestTypeWithSecondary(), value="test.ext"))
        inpsdict = WdlTranslator().build_inputs_file(w, merge_resources=False)
        self.assertEqual("test.ext", inpsdict.get("tst.wsec"))
        self.assertEqual("test.txt", inpsdict.get("tst.wsec_txt"))

    def test_inputs_generator_array_of_secondary_files(self):
        w = Workflow("tst")
        w._add_input(Input("wsec", Array(TestTypeWithSecondary()), value=["test.ext"]))
        inpsdict = WdlTranslator().build_inputs_file(w, merge_resources=False)
        self.assertListEqual(["test.ext"], inpsdict.get("tst.wsec"))
        self.assertListEqual(["test.txt"], inpsdict.get("tst.wsec_txt"))


class TestWdlSelectorsAndGenerators(unittest.TestCase):

    def test_input_selector_base_stringenv(self):
        input_sel = InputSelector("random")
        self.assertEqual("${random}", wdl.translate_input_selector(input_sel, None, string_environment=True))

    def test_input_selector_base_nostringenv(self):
        input_sel = InputSelector("random")
        self.assertEqual("random", wdl.translate_input_selector(input_sel, None, string_environment=False))

    def test_input_value_none_stringenv(self):
        self.assertEqual(None, wdl.get_input_value_from_potential_selector_or_generator(None, None, string_environment=True))

    def test_input_value_none_nostringenv(self):
        self.assertEqual(None, wdl.get_input_value_from_potential_selector_or_generator(None, None, string_environment=False))

    def test_input_value_string_stringenv(self):
        self.assertEqual(
            "TestString",
            wdl.get_input_value_from_potential_selector_or_generator("TestString", None, string_environment=True)
        )

    def test_input_value_string_nostringenv(self):
        self.assertEqual(
            '"TestString"',
            wdl.get_input_value_from_potential_selector_or_generator("TestString", None, string_environment=False)
        )

    def test_input_value_int_stringenv(self):
        self.assertEqual(
            42,
            wdl.get_input_value_from_potential_selector_or_generator(42, None, string_environment=True)
        )

    def test_input_value_int_nostringenv(self):
        self.assertEqual(
            42,
            wdl.get_input_value_from_potential_selector_or_generator(42, None, string_environment=False)
        )

    def test_input_value_filename_stringenv(self):
        import uuid
        fn = Filename(guid=str(uuid.uuid4()))
        self.assertEqual(
            fn.generated_filename(),
            wdl.get_input_value_from_potential_selector_or_generator(fn, None, string_environment=True)
        )

    def test_input_value_filename_nostringenv(self):
        import uuid
        fn = Filename(guid=str(uuid.uuid4()))
        self.assertEqual(
            '"%s"' % fn.generated_filename(),
            wdl.get_input_value_from_potential_selector_or_generator(fn, None, string_environment=False)
        )

    def test_input_value_wildcard(self):
        self.assertRaises(
            Exception,
            wdl.get_input_value_from_potential_selector_or_generator,
            value=WildcardSelector("*"),
            tool_id=None
        )

    def test_input_value_cpuselect_stringenv(self):
        inp = CpuSelector()
        self.assertEqual(
            "${runtime_cpu}",
            wdl.get_input_value_from_potential_selector_or_generator(inp, None, string_environment=True)
        )

    def test_input_value_cpuselect_nostringenv(self):
        inp = CpuSelector()
        self.assertEqual(
            "runtime_cpu",
            wdl.get_input_value_from_potential_selector_or_generator(inp, None, string_environment=False)
        )

    # def test_input_value_memselect_stringenv(self):
    #     inp = MemorySelector()
    #     self.assertEqual(
    #         "${floor(runtime_memory)}",
    #         wdl.get_input_value_from_potential_selector_or_generator(inp, string_environment=True)
    #     )
    #
    # def test_input_value_memselect_nostringenv(self):
    #     inp = MemorySelector()
    #     self.assertEqual(
    #         "floor(runtime_memory)",
    #         wdl.get_input_value_from_potential_selector_or_generator(inp, string_environment=False)
    #     )

    def test_input_value_wdl_callable(self):
        class CallableWdl:
            def wdl(self):
                return "unbelievable"

        self.assertEqual(
            "unbelievable",
            wdl.get_input_value_from_potential_selector_or_generator(CallableWdl(), None)
        )

    def test_input_value_wdl_noncallable(self):

        class NonCallableWdl:
            def __init__(self):
                self.wdl = None

        self.assertRaises(
            Exception,
            wdl.get_input_value_from_potential_selector_or_generator,
            value=NonCallableWdl(),
            tool_id=None
        )

    def test_string_formatter(self):
        b = StringFormatter("no format")
        res = wdl.get_input_value_from_potential_selector_or_generator(b, None)
        self.assertEqual("no format", res)

    def test_string_formatter_one_string_param(self):
        b = StringFormatter("there's {one} arg", one="a string")
        res = wdl.get_input_value_from_potential_selector_or_generator(b, None)
        self.assertEqual("there's a string arg", res)

    def test_string_formatter_one_input_selector_param(self):
        b = StringFormatter("an input {arg}", arg=InputSelector("random_input"))
        res = wdl.get_input_value_from_potential_selector_or_generator(b, None)
        self.assertEqual("an input ${random_input}", res)

    def test_string_formatter_two_param(self):
        # vardict input format
        b = StringFormatter("{tumorName}:{normalName}",
                            tumorName=InputSelector("tumorInputName"), normalName=InputSelector("normalInputName"))
        res = wdl.get_input_value_from_potential_selector_or_generator(b, None)
        self.assertEqual("${tumorInputName}:${normalInputName}", res)


class TestWdlGenerateInput(unittest.TestCase):

    def setUp(self):
        self.translator = wdl.WdlTranslator()

    def test_input_in_input_value_includetrue_nooptional_nodefault(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(), value="1", include_in_inputs_file_if_none=True))

        self.assertDictEqual({"test_input_in_inputfile.inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includetrue_nooptional_default(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(), value="1", default="2", include_in_inputs_file_if_none=True))

        self.assertDictEqual({"test_input_in_inputfile.inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includetrue_optional_nodefault(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(optional=True), value="1", include_in_inputs_file_if_none=True))

        self.assertDictEqual({"test_input_in_inputfile.inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includetrue_optional_default(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(optional=True), value="1", default="2", include_in_inputs_file_if_none=True))

        self.assertDictEqual({"test_input_in_inputfile.inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includefalse_nooptional_nodefault(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(), value="1", include_in_inputs_file_if_none=False))

        self.assertDictEqual({"test_input_in_inputfile.inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includefalse_nooptional_default(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(), value="1", default="2", include_in_inputs_file_if_none=False))

        self.assertDictEqual({"test_input_in_inputfile.inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includefalse_optional_nodefault(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(optional=True), value="1", include_in_inputs_file_if_none=False))

        self.assertDictEqual({"test_input_in_inputfile.inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includefalse_optional_default(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(
            Input("inpId", String(optional=True), value="1", default="2", include_in_inputs_file_if_none=False))

        self.assertDictEqual({"test_input_in_inputfile.inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includetrue_nooptional_nodefault(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(), include_in_inputs_file_if_none=True))

        self.assertDictEqual({"test_input_in_inputfile.inpId": None}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includetrue_nooptional_default(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(), default="2", include_in_inputs_file_if_none=True))

        self.assertDictEqual({"test_input_in_inputfile.inpId": None}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includetrue_optional_nodefault(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(optional=True), include_in_inputs_file_if_none=True))

        self.assertDictEqual({"test_input_in_inputfile.inpId": None}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includetrue_optional_default(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(optional=True), default="2", include_in_inputs_file_if_none=True))

        self.assertDictEqual({"test_input_in_inputfile.inpId": None}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includefalse_nooptional_nodefault(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(), include_in_inputs_file_if_none=False))

        # included because no value, no default, and not optional
        self.assertDictEqual({"test_input_in_inputfile.inpId": None}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includefalse_nooptional_default(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(), default="2", include_in_inputs_file_if_none=False))

        self.assertDictEqual({}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includefalse_optional_nodefault(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(optional=True), include_in_inputs_file_if_none=False))

        self.assertDictEqual({}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includefalse_optional_default(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(optional=True), default="2", include_in_inputs_file_if_none=False))

        self.assertDictEqual({}, self.translator.build_inputs_file(wf))

    def test_tool_output_with_input_selector(self):
        tool = TestToolWithSecondaryOutput()
        toolout = tool.outputs()[0]
        os = wdl.translate_output_node_with_glob(toolout, toolout.glob, tool)

        self.assertEqual("out", os[0].name)
        self.assertEqual('"${testtool}/out"', os[0].expression)

        self.assertEqual("out_txt", os[1].name)
        self.assertEqual('"${testtool}/out.txt"', os[1].expression)
