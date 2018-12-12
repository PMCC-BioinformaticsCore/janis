from Pipeline.utils.registry import Registry
from Pipeline.tool.tool import Tool

from typing import Type, List, Optional

__tools_registry = Registry[Type[Tool]]()


def register_tool(tool: Type[Tool]) -> bool:
    tool_id: str = tool.tool().lower()
    return __tools_registry.register(tool_id, tool)


def get_tool(tool_name: str) -> Optional[Type[Tool]]:
    return __tools_registry.get(tool_name.lower())


def get_tools() -> List[Type[Tool]]:
    return __tools_registry.objects()
