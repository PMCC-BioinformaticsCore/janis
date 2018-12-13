from Pipeline.utils.registry import Registry
from Pipeline.tool.commandtool import CommandTool

from typing import Type, List, Optional

__tools_registry = Registry[Type[CommandTool]]()


def register_tool(tool: Type[CommandTool]) -> bool:
    tool_id: str = tool.tool().lower()
    return __tools_registry.register(tool_id, tool)


def get_tool(tool_name: str) -> Optional[Type[CommandTool]]:
    return __tools_registry.get(tool_name.lower())


def get_tools() -> List[Type[CommandTool]]:
    return __tools_registry.objects()
