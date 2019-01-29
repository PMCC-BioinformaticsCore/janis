from janis.tool.commandtool import ToolInput
from janis.types.common_data_types import Array, File
from janis.unix.tools.unixtool import UnixTool


class Cat(UnixTool):

    @staticmethod
    def tool():
        return "cat"

    def friendly_name(self):
        return "Concatenate"

    @staticmethod
    def base_command():
        return "cat"

    def inputs(self):
        return [
            ToolInput("files", Array(File()))
        ]

    def outputs(self):
        return []
