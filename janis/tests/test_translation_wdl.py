import unittest
import wdlgen

from janis import ToolOutput, ToolInput, String, CommandTool, Stdout, InputSelector, Array, File


class TestWdl(unittest.TestCase):

    def test_optional_array(self):
        t = Array(File(), optional=True)
        wdl = t.wdl()
        self.assertIsInstance(wdl, wdlgen.WdlType)
        self.assertTrue(wdl.optional)
        self.assertEqual("Array[File]?", wdl.get_string())
