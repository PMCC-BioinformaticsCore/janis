import unittest
from typing import List

from janis.types import CpuSelector, MemorySelector

from janis import Workflow, ToolOutput, ToolInput, String, CommandTool, Stdout, InputSelector, Array, File, Filename, \
    WildcardSelector, Input, Output
import janis.translations.cwl as cwl


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


class TestCwlTypesConversion(unittest.TestCase):
    pass


class TestCwlMisc(unittest.TestCase):
    def test_str_tool(self):
        t = TestTool()
        self.assertEqual(t.translate("cwl"), cwl_testtool)


class TestCwlSelectorsAndGenerators(unittest.TestCase):

    def test_input_selector_base(self):
        input_sel = InputSelector("random")
        self.assertEqual("$(inputs.random)", cwl.translate_input_selector(input_sel))

    def test_input_selector_prefix(self):
        input_sel = InputSelector("random", prefix="&& ")
        self.assertEqual("&& $(inputs.random)", cwl.translate_input_selector(input_sel))

    def test_base_input_selector(self):
        input_sel = InputSelector("random", suffix=".cwl")
        self.assertEqual("$(inputs.random).cwl", cwl.translate_input_selector(input_sel))

    def test_input_value_none_stringenv(self):
        self.assertEqual(None, cwl.get_input_value_from_potential_selector_or_generator(None, "tool_id", string_environment=True))

    def test_input_value_none_nostringenv(self):
        self.assertEqual(None, cwl.get_input_value_from_potential_selector_or_generator(None, "tool_id", string_environment=False))

    def test_input_value_string_stringenv(self):
        self.assertEqual(
            "TestString",
            cwl.get_input_value_from_potential_selector_or_generator("TestString", "tool_id", string_environment=True)
        )

    def test_input_value_string_nostringenv(self):
        self.assertEqual(
            '"TestString"',
            cwl.get_input_value_from_potential_selector_or_generator("TestString", "tool_id", string_environment=False)
        )

    def test_input_value_int_stringenv(self):
        self.assertEqual(
            42,
            cwl.get_input_value_from_potential_selector_or_generator(42, "tool_id", string_environment=True)
        )

    def test_input_value_int_nostringenv(self):
        self.assertEqual(
            42,
            cwl.get_input_value_from_potential_selector_or_generator(42, "tool_id", string_environment=False)
        )

    def test_input_value_filename_stringenv(self):
        import uuid
        fn = Filename(guid=str(uuid.uuid4()))
        self.assertEqual(
            fn.generated_filename(),
            cwl.get_input_value_from_potential_selector_or_generator(fn, "tool_id", string_environment=True)
        )

    def test_input_value_filename_nostringenv(self):
        import uuid
        fn = Filename(guid=str(uuid.uuid4()))
        self.assertEqual(
            '"%s"' % fn.generated_filename(),
            cwl.get_input_value_from_potential_selector_or_generator(fn, "tool_id", string_environment=False)
        )

    def test_input_value_inpselect_stringenv(self):
        inp = InputSelector("threads")
        self.assertEqual(
            "$(inputs.threads)",
            cwl.get_input_value_from_potential_selector_or_generator(inp, "tool_id", string_environment=True)
        )

    def test_input_value_inpselect_nostringenv(self):
        inp = InputSelector("threads")
        self.assertEqual(
            "$(inputs.threads)",
            cwl.get_input_value_from_potential_selector_or_generator(inp, "tool_id", string_environment=False)
        )

    def test_input_value_wildcard(self):
        self.assertRaises(
            Exception,
            cwl.get_input_value_from_potential_selector_or_generator,
            value=WildcardSelector("*"),
            tool_id=None
        )

    def test_input_value_cpuselect_stringenv(self):
        inp = CpuSelector()
        self.assertEqual(
            "$(inputs.runtime_cpu)",
            cwl.get_input_value_from_potential_selector_or_generator(inp, "tool_id", string_environment=True)
        )

    def test_input_value_cpuselect_nostringenv(self):
        inp = CpuSelector()
        self.assertEqual(
            "$(inputs.runtime_cpu)",
            cwl.get_input_value_from_potential_selector_or_generator(inp, "tool_id", string_environment=False)
        )

    def test_input_value_memselect_stringenv(self):
        inp = MemorySelector()
        self.assertEqual(
            "$(Math.floor(inputs.runtime_memory))",
            cwl.get_input_value_from_potential_selector_or_generator(inp, "tool_id", string_environment=True)
        )

    def test_input_value_memselect_nostringenv(self):
        inp = MemorySelector()
        self.assertEqual(
            "$(Math.floor(inputs.runtime_memory))",
            cwl.get_input_value_from_potential_selector_or_generator(inp, "tool_id", string_environment=False)
        )

    def test_input_value_cwl_callable(self):
        class NonCallableCwl:
            def cwl(self):
                return "unbelievable"

        self.assertEqual(
            "unbelievable",
            cwl.get_input_value_from_potential_selector_or_generator(NonCallableCwl(), "tool_id")
        )

    def test_input_value_cwl_noncallable(self):

        class NonCallableCwl:
            def __init__(self):
                self.cwl = None

        self.assertRaises(
            Exception,
            cwl.get_input_value_from_potential_selector_or_generator,
            value=NonCallableCwl(),
            tool_id=None
        )


class TestCwlTranslateInput(unittest.TestCase):

    def test_translate_input(self):
        inp = Input(identifier="testIdentifier", data_type=String(), value="value",
                    label="myLabel", doc="docstring", default="defaultValue")
        tinp = cwl.translate_input(inp)

        self.assertEqual("testIdentifier", tinp.id)
        self.assertEqual("myLabel", tinp.label)
        self.assertEqual(None, tinp.secondaryFiles)
        self.assertEqual("docstring", tinp.doc)
        self.assertEqual(None, tinp.inputBinding)
        self.assertEqual("string", tinp.type)
        self.assertEqual("defaultValue", tinp.default)

    def test_secondary_file_translation(self):
        inp = Input(identifier="testIdentifier", data_type=TestTypeWithSecondary())
        tinp = cwl.translate_input(inp)

        self.assertEqual("File", tinp.type)
        self.assertListEqual([".txt"], tinp.secondaryFiles)


class TestCwlGenerateInput(unittest.TestCase):

    def setUp(self):
        self.translator = cwl.CwlTranslator()

    def test_input_in_input_value_includetrue_nooptional_nodefault(self):
        wf = Workflow("test_cwl_input_in_input_value_includetrue_nooptional_nodefault")
        wf.add_items(Input("inpId", String(), value="1", include_in_inputs_file_if_none=True))

        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includetrue_nooptional_default(self):
        wf = Workflow("test_cwl_input_in_input_value_includetrue_nooptional_default")
        wf.add_items(Input("inpId", String(), value="1", default="2", include_in_inputs_file_if_none=True))

        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includetrue_optional_nodefault(self):
        wf = Workflow("test_cwl_input_in_input_value_includetrue_optional_nodefault")
        wf.add_items(Input("inpId", String(optional=True), value="1", include_in_inputs_file_if_none=True))

        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includetrue_optional_default(self):
        wf = Workflow("test_cwl_input_in_input_value_includetrue_optional_default")
        wf.add_items(Input("inpId", String(optional=True), value="1", default="2", include_in_inputs_file_if_none=True))

        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includefalse_nooptional_nodefault(self):
        wf = Workflow("test_cwl_input_in_input_value_includefalse_nooptional_nodefault")
        wf.add_items(Input("inpId", String(), value="1", include_in_inputs_file_if_none=False))

        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includefalse_nooptional_default(self):
        wf = Workflow("test_cwl_input_in_input_value_includefalse_nooptional_default")
        wf.add_items(Input("inpId", String(), value="1", default="2", include_in_inputs_file_if_none=False))

        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includefalse_optional_nodefault(self):
        wf = Workflow("test_cwl_input_in_input_value_includefalse_optional_nodefault")
        wf.add_items(Input("inpId", String(optional=True), value="1", include_in_inputs_file_if_none=False))

        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_value_includefalse_optional_default(self):
        wf = Workflow("test_cwl_input_in_input_value_includefalse_optional_default")
        wf.add_items(
            Input("inpId", String(optional=True), value="1", default="2", include_in_inputs_file_if_none=False))

        self.assertDictEqual({"inpId": "1"}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includetrue_nooptional_nodefault(self):
        wf = Workflow("test_cwl_input_in_input_novalue_includetrue_nooptional_nodefault")
        wf.add_items(Input("inpId", String(), include_in_inputs_file_if_none=True))

        self.assertDictEqual({"inpId": None}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includetrue_nooptional_default(self):
        wf = Workflow("test_cwl_input_in_input_novalue_includetrue_nooptional_default")
        wf.add_items(Input("inpId", String(), default="2", include_in_inputs_file_if_none=True))

        self.assertDictEqual({"inpId": None}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includetrue_optional_nodefault(self):
        wf = Workflow("test_cwl_input_in_input_novalue_includetrue_optional_nodefault")
        wf.add_items(Input("inpId", String(optional=True), include_in_inputs_file_if_none=True))

        self.assertDictEqual({"inpId": None}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includetrue_optional_default(self):
        wf = Workflow("test_cwl_input_in_input_novalue_includetrue_optional_default")
        wf.add_items(Input("inpId", String(optional=True), default="2", include_in_inputs_file_if_none=True))

        self.assertDictEqual({"inpId": None}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includefalse_nooptional_nodefault(self):
        wf = Workflow("test_cwl_input_in_input_novalue_includefalse_nooptional_nodefault")
        wf.add_items(Input("inpId", String(), include_in_inputs_file_if_none=False))

        # included because no value, no default, and not optional
        self.assertDictEqual({"inpId": None}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includefalse_nooptional_default(self):
        wf = Workflow("test_cwl_input_in_input_novalue_includefalse_nooptional_default")
        wf.add_items(Input("inpId", String(), default="2", include_in_inputs_file_if_none=False))

        self.assertDictEqual({}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includefalse_optional_nodefault(self):
        wf = Workflow("test_cwl_input_in_input_novalue_includefalse_optional_nodefault")
        wf.add_items(Input("inpId", String(optional=True), include_in_inputs_file_if_none=False))

        self.assertDictEqual({}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includefalse_optional_default(self):
        wf = Workflow("test_cwl_input_in_input_novalue_includefalse_optional_default")
        wf.add_items(Input("inpId", String(optional=True), default="2", include_in_inputs_file_if_none=False))

        self.assertDictEqual({}, self.translator.build_inputs_file(wf))







cwl_testtool = """\
baseCommand: echo
class: CommandLineTool
cwlVersion: v1.0
id: testtranslation-tool
inputs:
- id: testtool
  label: testtool
  type: string
label: testtranslation-tool
outputs:
- id: std
  label: std
  type: stdout
requirements:
  DockerRequirement:
    dockerPull: ubuntu:latest
  InlineJavascriptRequirement: {}
"""
