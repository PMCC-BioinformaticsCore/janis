import unittest

from janis import Array, String


class TestTypes(unittest.TestCase):

    def test_array_of_strings(self):

        ar = Array(String())
        d = ar.cwl_type()
        self.assertEqual(d.get_dict(), {'type': 'array', 'items': 'string'})

    def test_array_of_array_of_strings(self):
        ar = Array(Array(String()))
        d = ar.cwl_type()
        self.assertEqual(d.get_dict(), {'type': 'array', 'items': {'type': 'array', 'items': 'string'}})
