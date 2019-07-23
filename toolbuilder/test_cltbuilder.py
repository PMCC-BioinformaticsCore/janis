import unittest

from toolbuilder.cltconvert import get_string_repr, convert_commandtool
from janis_core import InputSelector, ToolOutput, String, File
from janis_unix.data_types.csv import Csv
from janis_unix.tools.echo import Echo


class DataTypes(unittest.TestCase):
    def test_string_nooptional(self):
        d = String()
        self.assertEqual("String()", get_string_repr(d))

    def test_string_optional(self):
        d = String(optional=True)
        self.assertEqual("String(optional=True)", get_string_repr(d))

    def test_file_nooptional(self):
        d = File()
        self.assertEqual("File()", get_string_repr(d))

    def test_file_optional(self):
        d = File(optional=True)
        self.assertEqual("File(optional=True)", get_string_repr(d))

    def test_csv_nooptional(self):
        d = Csv(optional=False)
        self.assertEqual("Csv()", get_string_repr(d))

    def test_csv_optional(self):
        d = Csv(optional=True)
        self.assertEqual("Csv(optional=True)", get_string_repr(d))


class TestInputSelector(unittest.TestCase):
    def test_inputselector_1param(self):
        i = InputSelector("filename")
        self.assertEqual(
            'InputSelector(input_to_select="filename")', get_string_repr(i)
        )


class TestToolOutput(unittest.TestCase):
    def test_tooloutput(self):
        o = ToolOutput("tag", String())
        self.assertEqual(
            'ToolOutput(tag="tag", output_type=String())', get_string_repr(o)
        )


class TestEcho(unittest.TestCase):
    def test_echo(self):
        tool = Echo()
        print(convert_commandtool(tool))


class TestDateTime(unittest.TestCase):
    def test_datetime(self):
        from datetime import datetime

        dt = datetime(1985, 11, 26, 1, 20)
        self.assertEqual(
            'datetime.fromisoformat("1985-11-26T01:20:00")', get_string_repr(dt)
        )
