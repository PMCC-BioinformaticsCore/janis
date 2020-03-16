from typing import List

from janis_core import CodeTool
from .commandtool import prepare_commandtool_page


def prepare_code_tool_page(tool: CodeTool, versions: List[str]):
    return prepare_commandtool_page(tool, versions)
