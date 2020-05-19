import re
from typing import Tuple, List, Union
from enum import Enum

from .tool import tool_template, tool_version_template
from .gatk import gatk4_tool_template


class ToolTemplateType(Enum):
    base = "regular"
    gatk4 = "gatk4"


def generate_gatk4_tooltemplatebase(gatk_command, inputs, outputs, metadata):
    io_prefix = " " * 12

    return gatk4_tool_template.format(
        gatk_command=gatk_command,
        inputs=",\n".join(io_prefix + s for s in inputs),
        outputs=",\n".join(io_prefix + s for s in outputs),
        metadata=metadata,
    )


def generate_regular_tooltemplatebase(
    toolname: str,
    name: str,
    friendly_name: str,
    tool_provider: str,
    base_command: Union[str, List[str]],
    inputs: List[str],  # =",\n".join((io_prefix + get_string_repr(i)) for i in ins),
    outputs: List[str],  # =",\n".join((io_prefix + get_string_repr(o)) for o in outs),
    metadata: str,  # =get_string_repr(metadata),
):

    io_prefix = " " * 12

    return tool_template.format(
        toolname=toolname,
        name=name,
        friendly_name=friendly_name,
        tool_provider=tool_provider,
        base_command=base_command,
        inputs=",\n".join(io_prefix + s for s in inputs),
        outputs=",\n".join(io_prefix + s for s in outputs),
        metadata=metadata,
    )


def generate_tool_version_template(name: str, version: str, container: str):
    escapedversion = re.sub("[\W]", "_", str(version)) if version else str(None)

    return tool_version_template.format(
        name=name, container=container, version=version, escapedversion=escapedversion,
    )
