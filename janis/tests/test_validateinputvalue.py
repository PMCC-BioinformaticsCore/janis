import unittest
from janis import Input, String


class TestValidateInputvalue(unittest.TestCase):
    def test_validate_string_optional_allowoptional_value(self):
        v = Input("aa", String(optional=True), "value")
        self.assertEqual(True, v.validate_value(True))

    def test_validate_string_optional_allowoptional_novalue(self):
        v = Input("aa", String(optional=True))
        self.assertEqual(True, v.validate_value(True))

    def test_validate_string_optional_disallowoptional_value(self):
        v = Input("aa", String(optional=True), "value")
        self.assertEqual(True, v.validate_value(False))

    def test_validate_string_optional_disallowoptional_novalue(self):
        v = Input("aa", String(optional=True))
        self.assertEqual(True, v.validate_value(False))

    def test_validate_string_nooptional_allowoptional_value(self):
        v = Input("aa", String(), "value")
        self.assertEqual(True, v.validate_value(True))

    def test_validate_string_nooptional_allowoptional_novalue(self):
        v = Input("aa", String())
        self.assertEqual(True, v.validate_value(True))

    def test_validate_string_nooptional_disallowoptional_value(self):
        v = Input("aa", String(), "value")
        self.assertEqual(True, v.validate_value(False))

    def test_validate_string_nooptional_disallowoptional_novalue(self):
        v = Input("aa", String())
        self.assertEqual(False, v.validate_value(False))
