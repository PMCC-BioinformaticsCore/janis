import re
import inspect
from datetime import datetime
from typing import Union, List

from janis_core.tool.commandtool import ToolInput, ToolArgument, ToolOutput
from janis_core.types import (
    InputSelector,
    WildcardSelector,
    CpuSelector,
    MemorySelector,
)
from janis_core.types.data_types import DataType
from janis_core.utils.metadata import Metadata, ToolMetadata


tool_template = """
from datetime import datetime
from janis_core import CommandTool, ToolInput, ToolOutput, File, Boolean, String, Int, InputSelector, Filename, ToolMetadata

class {name}Base(CommandTool):

    def friendly_name(self) -> str:
        return "{friendly_name}"

    @staticmethod
    def tool_provider():
        return "{tool_provider}"

    @staticmethod
    def tool() -> str:
        return "{toolname}"

    @staticmethod
    def base_command():
        return {base_command}
    
    def inputs(self):
        return [
{inputs}
        ]
        
    def outputs(self):
        return [
{outputs}
        ]
        
    def metadata(self):
        return {metadata}
        
class {name}_{escapedversion}({name}Base):
    @staticmethod
    def version():
        return "{version}"
    
    @staticmethod
    def container():
        return "{container}"
"""


generic_convertible = [
    DataType,
    ToolInput,
    ToolOutput,
    ToolArgument,
    InputSelector,
    WildcardSelector,
    MemorySelector,
    CpuSelector,
    Metadata,
]


def get_string_repr(obj):

    if isinstance(obj, str):
        nlreplaced = obj.replace("\n", "\\n")
        return f'"{nlreplaced}"'
    elif isinstance(obj, datetime):
        return f'datetime.fromisoformat("{obj.isoformat()}")'
    elif any(isinstance(obj, T) for T in generic_convertible):
        return convert_generic_class(obj, type(obj).__name__)

    return str(obj)


def convert_generic_class(t, name, ignore_fields=None):
    options = []

    ignore_fields = set(ignore_fields if ignore_fields else ["self", "args", "kwargs"])
    params = inspect.signature(type(t).__init__).parameters
    # fields = fields_to_check if fields_to_check \
    #     else [f for f in dict(params).keys() if f not in ignore_fields]

    for fkey in params:
        if fkey in ignore_fields:
            continue
        opts = params[fkey]

        v = t.__getattribute__(fkey)
        if v is None and opts.default is None or v == opts.default:
            continue

        options.append(fkey + "=" + get_string_repr(v))

    return f"{name}({', '.join(options)})"


def convert_commandtool(commandtool):

    convert_command_tool_fragments(
        commandtool.id(),
        commandtool.base_command(),
        commandtool.friendly_name(),
        commandtool.tool_provider(),
        commandtool.inputs(),
        commandtool.outputs(),
        commandtool.metadata(),
        commandtool.version(),
        commandtool.docker(),
    )


def convert_command_tool_fragments(
    toolid: str,
    basecommand: Union[str, List[str]],
    friendly_name: str,
    toolprov: str,
    ins: [ToolInput],
    outs: [ToolOutput],
    metadata: ToolMetadata,
    version: str,
    container: str,
):
    io_prefix = "            "

    if not metadata.dateCreated:
        metadata.dateCreated = datetime.now()
    metadata.dateUpdated = datetime.now()

    if isinstance(toolid, list):
        toolid = "".join(s.title() for s in toolid)

    return tool_template.format(
        toolname=toolid,
        name=toolid.title(),
        friendly_name=friendly_name.title(),
        tool_provider=toolprov,
        base_command=f'"{basecommand}"'
        if isinstance(basecommand, str)
        else str(basecommand),
        inputs=",\n".join((io_prefix + get_string_repr(i)) for i in ins),
        outputs=",\n".join((io_prefix + get_string_repr(o)) for o in outs),
        metadata=get_string_repr(metadata),
        container=container,
        version=version,
        escapedversion=re.sub("[\W]", "_", str(version)) if version else str(None),
    )
