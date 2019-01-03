import unittest
from typing import Dict

from Pipeline.types.registry import register_type, get_type, get_types
from Pipeline.types.data_types import DataType, NativeType


class TestDataType(DataType):

    @staticmethod
    def name() -> str:
        return "TestDataType"

    @staticmethod
    def primitive() -> NativeType:
        return "test" # This will intentionally break in production

    @staticmethod
    def doc() -> str:
        return "This class is purely used to test how the register "


class TestRegistry(unittest.TestCase):

    def test_register(self):
        register_type(TestDataType)
        self.assertIsNotNone(get_type(TestDataType.name()))

    def test_case_indifference(self):
        register_type(TestDataType)
        inverted_case = "".join(a.upper() if a.islower() else a.lower() for a in TestDataType.name())
        self.assertIsNotNone(get_type(inverted_case))
