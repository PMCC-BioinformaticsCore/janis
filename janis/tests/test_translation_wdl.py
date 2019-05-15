import unittest
import wdlgen
from typing import List

from janis.types import CpuSelector, MemorySelector

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


class TestTypeWithSecondary(File):
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


class TestWdlSelectorsAndGenerators(unittest.TestCase):

    def test_input_selector_base_stringenv(self):
        input_sel = InputSelector("random")
        self.assertEqual("${random}", wdl.translate_input_selector(input_sel, string_environment=True))

    def test_input_selector_base_nostringenv(self):
        input_sel = InputSelector("random")
        self.assertEqual("random", wdl.translate_input_selector(input_sel, string_environment=False))

    def test_input_selector_prefix_stringenv(self):
        input_sel = InputSelector("random", prefix="&& ")
        self.assertEqual("&& ${random}", wdl.translate_input_selector(input_sel, string_environment=True))

    def test_input_selector_prefix_nostringenv(self):
        input_sel = InputSelector("random", prefix="&& ")
        self.assertEqual('"&& " + random', wdl.translate_input_selector(input_sel, string_environment=False))

    def test_base_input_selector_stringenv(self):
        input_sel = InputSelector("outputFilename", suffix=".wdl")
        self.assertEqual("${outputFilename}.wdl", wdl.translate_input_selector(input_sel, string_environment=True))

    def test_base_input_selector_nostringenv(self):
        input_sel = InputSelector("outputFilename", suffix=".wdl")
        self.assertEqual('outputFilename + ".wdl"', wdl.translate_input_selector(input_sel, string_environment=False))

    def test_input_value_none_stringenv(self):
        self.assertEqual(None, wdl.get_input_value_from_potential_selector_or_generator(None, "tool_id", string_environment=True))

    def test_input_value_none_nostringenv(self):
        self.assertEqual(None, wdl.get_input_value_from_potential_selector_or_generator(None, "tool_id", string_environment=False))

    def test_input_value_string_stringenv(self):
        self.assertEqual(
            "TestString",
            wdl.get_input_value_from_potential_selector_or_generator("TestString", "tool_id", string_environment=True)
        )

    def test_input_value_string_nostringenv(self):
        self.assertEqual(
            '"TestString"',
            wdl.get_input_value_from_potential_selector_or_generator("TestString", "tool_id", string_environment=False)
        )

    def test_input_value_int_stringenv(self):
        self.assertEqual(
            42,
            wdl.get_input_value_from_potential_selector_or_generator(42, "tool_id", string_environment=True)
        )

    def test_input_value_int_nostringenv(self):
        self.assertEqual(
            42,
            wdl.get_input_value_from_potential_selector_or_generator(42, "tool_id", string_environment=False)
        )

    def test_input_value_filename_stringenv(self):
        import uuid
        fn = Filename(guid=str(uuid.uuid4()))
        self.assertEqual(
            fn.generated_filename(),
            wdl.get_input_value_from_potential_selector_or_generator(fn, "tool_id", string_environment=True)
        )

    def test_input_value_filename_nostringenv(self):
        import uuid
        fn = Filename(guid=str(uuid.uuid4()))
        self.assertEqual(
            '"%s"' % fn.generated_filename(),
            wdl.get_input_value_from_potential_selector_or_generator(fn, "tool_id", string_environment=False)
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
            wdl.get_input_value_from_potential_selector_or_generator(inp, "tool_id", string_environment=True)
        )

    def test_input_value_cpuselect_nostringenv(self):
        inp = CpuSelector()
        self.assertEqual(
            "runtime_cpu",
            wdl.get_input_value_from_potential_selector_or_generator(inp, "tool_id", string_environment=False)
        )

    # def test_input_value_memselect_stringenv(self):
    #     inp = MemorySelector()
    #     self.assertEqual(
    #         "${floor(runtime_memory)}",
    #         wdl.get_input_value_from_potential_selector_or_generator(inp, "tool_id", string_environment=True)
    #     )
    #
    # def test_input_value_memselect_nostringenv(self):
    #     inp = MemorySelector()
    #     self.assertEqual(
    #         "floor(runtime_memory)",
    #         wdl.get_input_value_from_potential_selector_or_generator(inp, "tool_id", string_environment=False)
    #     )

    def test_input_value_wdl_callable(self):
        class CallableWdl:
            def wdl(self):
                return "unbelievable"

        self.assertEqual(
            "unbelievable",
            wdl.get_input_value_from_potential_selector_or_generator(CallableWdl(), "tool_id")
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


class TestWdlTranslateInputs(unittest.TestCase):

    def test_input_in_input_value_includetrue(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(), value="1", include_in_inputs_file_if_none=True))
        translator = wdl.WdlTranslator()

        self.assertDictEqual({"test_input_in_inputfile.inpId": "1"}, translator.build_inputs_file(wf))

    def test_input_in_input_value_includefalse(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(), value="1", include_in_inputs_file_if_none=False))
        translator = wdl.WdlTranslator()

        self.assertDictEqual({"test_input_in_inputfile.inpId": "1"}, translator.build_inputs_file(wf))

    def test_input_in_input_with_novalue_includetrue(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(), value=None, include_in_inputs_file_if_none=True))
        translator = wdl.WdlTranslator()

        self.assertDictEqual({"test_input_in_inputfile.inpId": None}, translator.build_inputs_file(wf))

    def test_input_in_input_with_novalue_includefalse(self):
        wf = Workflow("test_input_in_inputfile")
        wf.add_items(Input("inpId", String(), value=None, include_in_inputs_file_if_none=False))
        translator = wdl.WdlTranslator()

        self.assertDictEqual({}, translator.build_inputs_file(wf))
