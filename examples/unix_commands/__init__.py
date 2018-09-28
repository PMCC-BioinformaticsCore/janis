from pipeline_definition.types.type_registry import register_input_factory

from examples.unix_commands.data_types.tar_file import TarFileFactory
from examples.unix_commands.data_types.generic_file import GenericFileFactory

register_input_factory(TarFileFactory())
register_input_factory(GenericFileFactory())
