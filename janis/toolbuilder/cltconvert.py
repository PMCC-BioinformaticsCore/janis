import inspect
from datetime import datetime

from janis.tool.tool import ToolInput, ToolArgument, ToolOutput
from janis.types import InputSelector, WildcardSelector, CpuSelector, MemorySelector
from janis.types.data_types import DataType


tool_template = """
from datetime import datetime
from janis import CommandTool, ToolInput, ToolOutput, File, Boolean, String, Int, InputSelector, Filename

class {name}Base(CommandTool):

    def friendly_name(self) -> str:
        return "{friendly_name}"

    @staticmethod
    def tool_provider():
        return "{tool_provider}"

    @staticmethod
    def tool() -> str:
        return "{name}"

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
        return ToolMetadata(
            creator=None, 
            maintainer=None, maintainer_email=None,
            date_created=datetime({cy}, {cm}, {cd}), date_updated=datetime({cy}, {cm}, {cd}),
            institution=None, doi=None,
            citation=None,
            keywords=["{name}"],
            documentation_url="{url}",
            documentation=\"""{doc}""\")
"""


generic_convertible = [DataType, ToolInput, ToolOutput, ToolArgument, InputSelector, WildcardSelector, MemorySelector, CpuSelector]


def get_string_repr(obj):

    if isinstance(obj, str):
        return f'"{obj}"'
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
        if fkey in ignore_fields: continue
        opts = params[fkey]

        v = t.__getattribute__(fkey)
        if v is None and opts.default is None or v == opts.default: continue

        options.append(fkey + "=" + get_string_repr(v))

    return f"{name}({', '.join(options)})"


def convert_commandtool(commandtool):

    ins = [get_string_repr(ti) for ti in commandtool.inputs()]
    outs = [get_string_repr(to) for to in commandtool.outputs()]

    bc = commandtool.base_command()
    io_prefix = "            "

    return tool_template.format(
        name=commandtool.id(),
        friendly_name=commandtool.friendly_name(),
        base_command=f'"{bc}"' if isinstance(bc, str) else str(bc),
        inputs=",\n".join((io_prefix + i) for i in ins),
        outputs=",\n".join((io_prefix + o) for o in outs),
        cy=datetime.now().year,
        cm=datetime.now().month,
        cd=datetime.now().day,
        doc=commandtool.metadata().documentation,
        url=commandtool.metadata().documentationUrl
    )


# def convert_toolargument(toolargument) -> str:
#     return convert_generic_class(toolargument, "ToolArgument")


# def convert_toolinput(toolinput) -> str:
#     options = []
#
#     options_to_check = ["prefix", "position", "shell_quote", "separate_value_from_prefix", "default",
#                         "nest_input_binding_on_array", "separator", "localise_file", "doc"]
#
#     for o in options_to_check:
#         v = toolinput.__getattribute__(o)
#         if isinstance(v, str):
#             v = f'"{v}"'
#         elif isinstance(v, InputSelector):
#             v = "InputSelec"
#         options.append(o + "=" + v)
#
#     return toolinout_template.format(
#         tooltype="Input",
#         id=toolinput.id(),
#         type=type(toolinput.input_type).__name__,
#         optionality="optional=True" if toolinput.input_type.optional else "",
#         options="".join(", " + o for o in options)
#     )
