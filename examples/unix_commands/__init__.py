from pipeline_definition.types.type_registry import register_tool, register_type

from examples.unix_commands.data_types.tar_file import TarFile

from examples.unix_commands.steps.compile import Compile
from examples.unix_commands.steps.untar import Untar
from examples.unix_commands.steps.tar import Tar


register_type(TarFile)

register_tool(Tar)
register_tool(Untar)
register_tool(Compile)

