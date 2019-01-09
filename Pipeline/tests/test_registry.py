import unittest
from typing import Dict, List

from Pipeline import CommandTool, ToolOutput, ToolInput, register_tool, get_tool
from Pipeline.types.registry import register_type, get_type, get_types
from Pipeline.types.data_types import DataType, NativeType
from Pipeline.utils.registry import TaggedRegistry


class TestDataType(DataType):

    @staticmethod
    def name() -> str:
        return "TestDataType"

    @staticmethod
    def primitive() -> NativeType:
        return "test"  # This will intentionally break in production

    @staticmethod
    def doc() -> str:
        return "This class is purely used to test how the register"


class TestTool(CommandTool):

    @staticmethod
    def tool():
        return "TestRegistryTool"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self) -> List[ToolInput]:
        return []

    def outputs(self) -> List[ToolOutput]:
        return []


class TestTaggedTool(TestTool):

    @staticmethod
    def version():
        return "1.0"


class TestTaggedToolHigher(TestTool):

    @staticmethod
    def version():
        return "1.1"


class TestRegistry(unittest.TestCase):

    def test_register(self):
        register_type(TestDataType)
        self.assertIsNotNone(get_type(TestDataType.name()))

    def test_case_indifference(self):
        register_type(TestDataType)
        inverted_case = "".join(a.upper() if a.islower() else a.lower() for a in TestDataType.name())
        self.assertIsNotNone(get_type(inverted_case))


class TestTaggedRegistry(unittest.TestCase):
    def test_register_nontagged(self):
        r = TaggedRegistry("latest")
        tool = TestTool()

        self.assertTrue(r.register(tool.tool(), tool.version(), tool))
        fetched = r.get(tool.tool(), None)
        self.assertIsNotNone(fetched)

    def test_register_tagged(self):
        r = TaggedRegistry("latest")
        tool = TestTaggedTool()

        self.assertTrue(r.register(tool.tool(), tool.version(), tool))
        fetched = r.get(tool.tool(), "1.0")
        self.assertIsNotNone(fetched)

    def test_register_tagged_multiple(self):
        r = TaggedRegistry("latest")

        tool_id = TestTool.tool()
        tool1 = TestTaggedTool()
        tool2 = TestTaggedToolHigher()

        self.assertTrue(r.register(tool_id, tool1.version(), tool1))
        self.assertTrue(r.register(tool_id, tool2.version(), tool2))
        self.assertEqual(tool1, r.get(tool_id, "1.0"))
        self.assertEqual(tool2, r.get(tool_id, "1.1"))

    def test_register_latest(self):
        r = TaggedRegistry("latest")
        tool = TestTool()

        self.assertTrue(r.register(tool.tool(), tool.version(), tool))
        fetched = r.get(tool.tool(), "latest")
        self.assertIsNotNone(fetched)

    def test_higher_version_override(self):
        r = TaggedRegistry("latest")
        tool_id = TestTool.tool()
        tool1 = TestTaggedTool()
        tool2 = TestTaggedToolHigher()

        self.assertTrue(r.register(tool_id, tool1.version(), tool1))
        self.assertEqual(tool1, r.get(tool_id, "latest"))
        self.assertTrue(r.register(tool_id, tool2.version(), tool2))
        self.assertEqual(tool2, r.get(tool_id, "latest"))

    def test_nontagged_vs_higher_override(self):
        r = TaggedRegistry("latest")
        tool_id = TestTool.tool()
        tool1 = TestTool()
        tool2 = TestTaggedToolHigher()

        self.assertTrue(r.register(tool_id, tool1.version(), tool1))
        self.assertEqual(tool1, r.get(tool_id, "latest"))
        self.assertTrue(r.register(tool_id, tool2.version(), tool2))
        self.assertEqual(tool2, r.get(tool_id, "latest"))

    def test_get_no_tag(self):
        r = TaggedRegistry("latest")

        tool = TestTool()
        self.assertTrue(r.register(tool.tool(), tool.version(), tool))
        self.assertTrue(tool, r.get(tool.tool(), None))
