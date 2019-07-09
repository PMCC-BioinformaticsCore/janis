import unittest
from os import getcwd
import os.path

from janis import String, Input
from janis.translations.exportpath import ExportPathKeywords
from janis.translations.translationbase import TranslatorBase


class TestExportPath(unittest.TestCase):
    """
    These tests are for components that apply to every single translation.
    Perhaps the export path function, or anything in translationbase.py.
    """

    def test_user_path(self):
        self.assertEqual(
            os.path.expanduser("~"), ExportPathKeywords.resolve("~", None, None)
        )

    def test_workflow_spec(self):
        self.assertEqual(
            "my/path/to/cwl",
            ExportPathKeywords.resolve("my/path/to/{language}", "cwl", None),
        )

    def test_workflow_name(self):
        self.assertEqual(
            "my/workflow_name/path",
            ExportPathKeywords.resolve("my/{name}/path", None, "workflow_name"),
        )

    def test_multi_replace(self):
        self.assertEqual(
            "test_multi_replace/test_multi_replace/test_multi_replace",
            ExportPathKeywords.resolve(
                "{name}/{name}/{name}", None, "test_multi_replace"
            ),
        )

    def test_combo_replace(self):
        self.assertEqual(
            os.path.expanduser("~") + "/Desktop/workflowname/wdl/",
            ExportPathKeywords.resolve(
                "~/Desktop/{name}/{language}/", "wdl", "workflowname"
            ),
        )

    def test_replace_only_cwd_dot(self):
        self.assertEqual(getcwd(), ExportPathKeywords.resolve(".", None, None))

    def test_replace_none(self):
        self.assertEqual(getcwd(), ExportPathKeywords.resolve(None, None, None))

    def test_replace_empty(self):
        self.assertEqual(getcwd(), ExportPathKeywords.resolve("", None, None))

    def test_replace_falsey(self):
        self.assertEqual(getcwd(), ExportPathKeywords.resolve(False, None, None))

    def test_replace_cwd_in_scope(self):
        self.assertEqual(
            os.path.join(getcwd(), "test"),
            ExportPathKeywords.resolve("./test", None, None),
        )

    def test_replace_dontreplace(self):
        path = ".myfile/starting/with/dot"
        self.assertEqual(path, ExportPathKeywords.resolve(path, None, None))

    def test_random_dot(self):
        path = "/mypath/./starting/with/dot"
        self.assertEqual(path, ExportPathKeywords.resolve(path, None, None))

    def test_no_spec_except(self):
        self.assertRaises(
            Exception,
            ExportPathKeywords.resolve,
            path="{language}",
            workflow_spec=None,
            workflow_name="name",
        )

    def test_no_name_except(self):
        self.assertRaises(
            Exception,
            ExportPathKeywords.resolve,
            path="{name}",
            workflow_spec="spec",
            workflow_name=None,
        )


class TestCanSkipIncludeTrueImport(unittest.TestCase):
    def test_input_in_input_value_includetrue_nooptional_nodefault(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input("inpId", String(), value="1", include_in_inputs_file_if_none=True)
        )
        self.assertEqual(False, can_skip)

    def test_input_in_input_value_includetrue_nooptional_default(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input(
                "inpId",
                String(),
                value="1",
                default="2",
                include_in_inputs_file_if_none=True,
            )
        )
        self.assertEqual(False, can_skip)

    def test_input_in_input_value_includetrue_optional_nodefault(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input(
                "inpId",
                String(optional=True),
                value="1",
                include_in_inputs_file_if_none=True,
            )
        )
        self.assertEqual(False, can_skip)

    def test_input_in_input_value_includetrue_optional_default(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input(
                "inpId",
                String(optional=True),
                value="1",
                default="2",
                include_in_inputs_file_if_none=True,
            )
        )
        self.assertEqual(False, can_skip)

    def test_input_in_input_value_includefalse_nooptional_nodefault(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input("inpId", String(), value="1", include_in_inputs_file_if_none=False)
        )
        self.assertEqual(False, can_skip)

    def test_input_in_input_value_includefalse_nooptional_default(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input(
                "inpId",
                String(),
                value="1",
                default="2",
                include_in_inputs_file_if_none=False,
            )
        )
        self.assertEqual(False, can_skip)

    def test_input_in_input_value_includefalse_optional_nodefault(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input(
                "inpId",
                String(optional=True),
                value="1",
                include_in_inputs_file_if_none=False,
            )
        )
        self.assertEqual(False, can_skip)

    def test_input_in_input_value_includefalse_optional_default(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input(
                "inpId",
                String(optional=True),
                value="1",
                default="2",
                include_in_inputs_file_if_none=False,
            )
        )
        self.assertEqual(False, can_skip)

    def test_input_in_input_novalue_includetrue_nooptional_nodefault(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input("inpId", String(), include_in_inputs_file_if_none=True)
        )
        self.assertEqual(False, can_skip)

    def test_input_in_input_novalue_includetrue_nooptional_default(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input("inpId", String(), default="2", include_in_inputs_file_if_none=True)
        )
        self.assertEqual(False, can_skip)

    def test_input_in_input_novalue_includetrue_optional_nodefault(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input("inpId", String(optional=True), include_in_inputs_file_if_none=True)
        )
        self.assertEqual(False, can_skip)

    def test_input_in_input_novalue_includetrue_optional_default(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input(
                "inpId",
                String(optional=True),
                default="2",
                include_in_inputs_file_if_none=True,
            )
        )
        self.assertEqual(False, can_skip)

    def test_input_in_input_novalue_includefalse_nooptional_nodefault(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input("inpId", String(), include_in_inputs_file_if_none=False)
        )
        self.assertEqual(False, can_skip)

    def test_input_in_input_novalue_includefalse_nooptional_default(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input("inpId", String(), default="2", include_in_inputs_file_if_none=False)
        )
        # new interpretation: can't skip because default
        self.assertEqual(False, can_skip)
        # self.assertEqual(True, can_skip)

    def test_input_in_input_novalue_includefalse_optional_nodefault(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input("inpId", String(optional=True), include_in_inputs_file_if_none=False)
        )
        self.assertEqual(True, can_skip)

    def test_input_in_input_novalue_includefalse_optional_default(self):
        can_skip = TranslatorBase.inp_can_be_skipped(
            Input(
                "inpId",
                String(optional=True),
                default="2",
                include_in_inputs_file_if_none=False,
            )
        )
        # new interpretation: can't skip because default
        self.assertEqual(False, can_skip)
        # self.assertEqual(True, can_skip)
