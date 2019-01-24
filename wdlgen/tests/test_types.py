import unittest

import wdlgen.types as types


class TestTypes(unittest.TestCase):

    def test_parse_string(self):
        t = types.WdlType.parse_type("String")
        self.assertIsInstance(t, types.WdlType)
        self.assertEqual(t._type._type, "String")
        self.assertFalse(t._optional)

    def test_parse_optional_string(self):
        t = types.WdlType.parse_type("String?")
        self.assertIsInstance(t, types.WdlType)
        self.assertEqual(t._type._type, "String")
        self.assertTrue(t._optional)

    def test_parse_primitives(self):
        results = [types.WdlType.parse_type(t) for t in types.PrimitiveType.types]
        failed = [r for r in results if r is None]
        self.assertEqual(len(failed), 0)

    def test_parse_array(self):
        t = types.WdlType.parse_type("Array[String]")

        tt = t._type
        self.assertIsInstance(tt, types.ArrayType)
        self.assertEqual(tt._subtype._type._type, "String")

    def test_parse_nested_array(self):
        t = types.WdlType.parse_type("Array[Array[Int]]")

        tt = t._type
        self.assertIsInstance(tt, types.ArrayType)
        ttt = tt._subtype._type
        self.assertIsInstance(ttt, types.ArrayType)
        self.assertEqual(ttt._subtype._type._type, "Int")

    def test_parse_qualified_array(self):
        t = types.WdlType.parse_type("Array[File]+")

        tt = t._type
        self.assertIsInstance(tt, types.ArrayType)
        self.assertTrue(tt._requires_multiple)


