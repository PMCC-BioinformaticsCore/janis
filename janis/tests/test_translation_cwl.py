import unittest
import cwlgen
from typing import List

from janis.translations import CwlTranslator

from janis.types import CpuSelector, MemorySelector

from janis import Workflow, ToolOutput, ToolInput, String, CommandTool, Stdout, InputSelector, Array, File, Filename, \
    WildcardSelector, Input, Output, StringFormatter
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


class TestToolWithSecondaryOutput(TestTool):
    def outputs(self):
        return [
            ToolOutput("out", TestTypeWithSecondary(), glob=InputSelector("testtool") + "/out")
        ]


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


class TestCwlTranslatorOverrides(unittest.TestCase):

    def setUp(self):
        self.translator = CwlTranslator()

    def test_stringify_workflow(self):
        cwlobj = cwlgen.Workflow("wid")
        self.assertEqual(
            "class: Workflow\ncwlVersion: v1.0\nid: wid\ninputs: {}\noutputs: {}\nsteps: {}\n",
            self.translator.stringify_translated_workflow(cwlobj)
        )

    def test_stringify_tool(self):
        cwlobj = cwlgen.CommandLineTool("tid")
        self.assertEqual(
            "class: CommandLineTool\ncwlVersion: v1.0\nid: tid\n",
            self.translator.stringify_translated_tool(cwlobj)
        )

    def test_stringify_inputs(self):
        d = {"inp1": 1}
        self.assertEqual(
            "inp1: 1\n",
            self.translator.stringify_translated_inputs(d)
        )

    def test_workflow_filename(self):
        w = Workflow("wid")
        self.assertEqual("wid.cwl", self.translator.workflow_filename(w))

    def test_tools_filename(self):
        self.assertEqual("TestTranslation-tool.cwl", self.translator.tool_filename(TestTool()))

    def test_inputs_filename(self):
        w = Workflow("wid")
        self.assertEqual("wid-inp.yml", self.translator.inputs_filename(w))

    def test_resources_filename(self):
        w = Workflow("wid")
        self.assertEqual("wid-resources.yml", self.translator.resources_filename(w))


class TestCwlArraySeparators(unittest.TestCase):
    # Based on https://www.commonwl.org/user_guide/09-array-inputs/index.html

    def test_regular_input_bindingin(self):
        t = ToolInput("filesA", Array(String()), prefix="-A", position=1)
        cwltoolinput = cwl.translate_tool_input(t)
        self.assertDictEqual({
            'id': 'filesA',
            'label': 'filesA',
            'type': {'items': 'string', 'type': 'array'},
            'inputBinding': {
                'prefix': '-A',
                'position': 1
            }
        }, cwltoolinput.get_dict())

    def test_nested_input_binding(self):
        t = ToolInput("filesB", Array(String()), prefix="-B=", separate_value_from_prefix=False,
                      position=2, prefix_applies_to_all_elements=True)
        cwltoolinput = cwl.translate_tool_input(t)
        self.assertDictEqual({
            'id': 'filesB',
            'label': 'filesB',
            'type': {
                'items': 'string',
                'type': 'array',
                'inputBinding': {
                    'prefix': '-B=',
                    'separate': False
                }
            },
            'inputBinding': {
                'position': 2
            }
        }, cwltoolinput.get_dict())

    def test_separated_input_bindingin(self):
        t = ToolInput("filesC", Array(String()), prefix="-C=", separate_value_from_prefix=False,
                      position=4, separator=",")
        cwltoolinput = cwl.translate_tool_input(t)
        self.assertDictEqual({
            'id': 'filesC',
            'label': 'filesC',
            'type': {
                'items': 'string',
                'type': 'array'
            },
            'inputBinding': {
                'prefix': '-C=',
                'itemSeparator': ',',
                'separate': False,
                'position': 4,
            }
        }, cwltoolinput.get_dict())

    def test_optional_array_prefixes(self):
        t = ToolInput("filesD", Array(String(), optional=True), prefix="-D", prefix_applies_to_all_elements=True)
        cwltoolinput = cwl.translate_tool_input(t)

        self.assertDictEqual({
            'id': 'filesD',
            'label': 'filesD',
            'type': [{
                'inputBinding': { 'prefix': '-D' },
                'items': 'string',
                'type': 'array'
            }, 'null'
            ]
        }, cwltoolinput.get_dict())


class TestCwlSelectorsAndGenerators(unittest.TestCase):

    def test_input_selector_base(self):
        input_sel = InputSelector("random")
        self.assertEqual("$(inputs.random)", cwl.translate_input_selector(input_sel, code_environment=False))\

    def test_input_selector_base_codeenv(self):
        input_sel = InputSelector("random")
        self.assertEqual("inputs.random", cwl.translate_input_selector(input_sel, code_environment=True))

    def test_input_value_none_codeenv(self):
        self.assertEqual(None, cwl.get_input_value_from_potential_selector_or_generator(None, code_environment=True))

    def test_input_value_none_nocodeenv(self):
        self.assertEqual(None, cwl.get_input_value_from_potential_selector_or_generator(None, code_environment=False))

    def test_input_value_string_codeenv(self):
        self.assertEqual(
            '"TestString"',
            cwl.get_input_value_from_potential_selector_or_generator("TestString", code_environment=True)
        )

    def test_input_value_string_nocodeenv(self):
        self.assertEqual(
            'TestString',
            cwl.get_input_value_from_potential_selector_or_generator("TestString", code_environment=False)
        )

    def test_input_value_int_codeenv(self):
        self.assertEqual(
            42,
            cwl.get_input_value_from_potential_selector_or_generator(42, code_environment=True)
        )

    def test_input_value_int_nocodeenv(self):
        self.assertEqual(
            42,
            cwl.get_input_value_from_potential_selector_or_generator(42, code_environment=False)
        )

    # def test_input_value_filename_codeenv(self):
    #     import uuid
    #     fn = Filename(guid=str(uuid.uuid4()))
    #     self.assertEqual(
    #         '"generated-" + Math.random().toString(16).substring(2, 8) + ""',
    #         cwl.get_input_value_from_potential_selector_or_generator(fn, code_environment=True)
    #     )
    #
    # def test_input_value_filename_nocodeenv(self):
    #     import uuid
    #     fn = Filename(guid=str(uuid.uuid4()))
    #     self.assertEqual(
    #         '$("generated-" + Math.random().toString(16).substring(2, 8) + "")',
    #         cwl.get_input_value_from_potential_selector_or_generator(fn, code_environment=False)
    #     )

    def test_input_value_inpselect_codeenv(self):
        inp = InputSelector("threads")
        self.assertEqual(
            "inputs.threads",
            cwl.get_input_value_from_potential_selector_or_generator(inp, code_environment=True)
        )

    def test_input_value_inpselect_nocodeenv(self):
        inp = InputSelector("threads")
        self.assertEqual(
            "$(inputs.threads)",
            cwl.get_input_value_from_potential_selector_or_generator(inp, code_environment=False)
        )

    def test_input_value_wildcard(self):
        self.assertRaises(
            Exception,
            cwl.get_input_value_from_potential_selector_or_generator,
            value=WildcardSelector("*")
        )

    def test_input_value_cpuselect_codeenv(self):
        inp = CpuSelector()
        self.assertEqual(
            "$(inputs.runtime_cpu)",
            cwl.get_input_value_from_potential_selector_or_generator(inp, code_environment=True)
        )

    def test_input_value_cpuselect_nocodeenv(self):
        inp = CpuSelector()
        self.assertEqual(
            "$(inputs.runtime_cpu)",
            cwl.get_input_value_from_potential_selector_or_generator(inp, code_environment=False)
        )

    def test_input_value_memselect_codeenv(self):
        inp = MemorySelector()
        self.assertEqual(
            "$(Math.floor(inputs.runtime_memory))",
            cwl.get_input_value_from_potential_selector_or_generator(inp, code_environment=True)
        )

    def test_input_value_memselect_nocodeenv(self):
        inp = MemorySelector()
        self.assertEqual(
            "$(Math.floor(inputs.runtime_memory))",
            cwl.get_input_value_from_potential_selector_or_generator(inp, code_environment=False)
        )

    def test_input_value_cwl_callable(self):
        class NonCallableCwl:
            def cwl(self):
                return "unbelievable"

        self.assertEqual(
            "unbelievable",
            cwl.get_input_value_from_potential_selector_or_generator(NonCallableCwl())
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

    def test_string_formatter(self):
        b = StringFormatter("no format")
        res = cwl.get_input_value_from_potential_selector_or_generator(b)
        self.assertEqual("no format", res)

    def test_string_formatter_one_string_param(self):
        b = StringFormatter("there's {one} arg", one="a string")
        res = cwl.get_input_value_from_potential_selector_or_generator(b)
        self.assertEqual('$("there\'s {one} arg".replace(/\{one\}/g, "a string"))', res)

    def test_string_formatter_one_input_selector_param(self):
        b = StringFormatter("an input {arg}", arg=InputSelector("random_input"))
        res = cwl.get_input_value_from_potential_selector_or_generator(b, code_environment=False)
        self.assertEqual('$("an input {arg}".replace(/\{arg\}/g, inputs.random_input))', res)

    def test_string_formatter_two_param(self):
        # vardict input format
        b = StringFormatter("{tumorName}:{normalName}",
                            tumorName=InputSelector("tumorInputName"), normalName=InputSelector("normalInputName"))
        res = cwl.get_input_value_from_potential_selector_or_generator(b)
        self.assertEqual('$("{tumorName}:{normalName}".replace(/\{tumorName\}/g, inputs.tumorInputName).replace(/\{normalName\}/g, inputs.normalInputName))', res)


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

        self.assertDictEqual({"inpId": None}, self.translator.build_inputs_file(wf))
        # self.assertDictEqual({}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includefalse_optional_nodefault(self):
        wf = Workflow("test_cwl_input_in_input_novalue_includefalse_optional_nodefault")
        wf.add_items(Input("inpId", String(optional=True), include_in_inputs_file_if_none=False))

        self.assertDictEqual({'inpId': None}, self.translator.build_inputs_file(wf))
        # self.assertDictEqual({}, self.translator.build_inputs_file(wf))

    def test_input_in_input_novalue_includefalse_optional_default(self):
        wf = Workflow("test_cwl_input_in_input_novalue_includefalse_optional_default")
        wf.add_items(Input("inpId", String(optional=True), default="2", include_in_inputs_file_if_none=False))

        # self.assertDictEqual({}, self.translator.build_inputs_file(wf))
        self.assertDictEqual({'inpId': None}, self.translator.build_inputs_file(wf))


cwl_testtool = """\
baseCommand: echo
class: CommandLineTool
cwlVersion: v1.0
id: TestTranslation-tool
inputs:
- id: testtool
  label: testtool
  type: string
label: TestTranslation-tool
outputs:
- id: std
  label: std
  type: stdout
requirements:
  DockerRequirement:
    dockerPull: ubuntu:latest
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
"""
