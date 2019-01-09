from typing import Type, List, Optional

from Pipeline.tool.commandtool import CommandTool
from Pipeline.tool.tool import Tool
from Pipeline.utils.registry import TaggedRegistry

__tools_registry = TaggedRegistry[Type[CommandTool]](default_tag="latest")


def register_tool(tool: Type[Tool], tool_id: Optional[str]=None) -> bool:
    if tool_id is None:
        tool_id = tool.tool().lower()

    version = None
    if callable(getattr(tool, "version", None)):
        version = tool.version()
    return __tools_registry.register(tool_id, version, tool)


def get_tool(tool_name: str, version: Optional[str]) -> Optional[Type[CommandTool]]:
    return __tools_registry.get(tool_name.lower(), version)


def get_tools() -> List[Type[CommandTool]]:
    return __tools_registry.objects()
