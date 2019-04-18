import unittest
from os.path import expanduser

from janis.translations.exportpath import ExportPathKeywords


class TestExportPath(unittest.TestCase):

    def test_user_path(self):
        from os.path import expanduser

        self.assertEqual(expanduser("~"), ExportPathKeywords.resolve("~", None, None))

    def test_workflow_spec(self):
        self.assertEqual("my/path/to/cwl", ExportPathKeywords.resolve("my/path/to/{language}", "cwl", None))

    def test_workflow_name(self):
        self.assertEqual("my/workflow_name/path", ExportPathKeywords.resolve("my/{name}/path", None, "workflow_name"))

    def test_multi_replace(self):
        self.assertEqual(
            "test_multi_replace/test_multi_replace/test_multi_replace",
            ExportPathKeywords.resolve("{name}/{name}/{name}", None, "test_multi_replace")
        )

    def test_combo_replace(self):
        self.assertEqual(
            expanduser('~') + "/Desktop/workflowname/wdl/",
            ExportPathKeywords.resolve("~/Desktop/{name}/{language}/", "wdl", "workflowname")
        )

