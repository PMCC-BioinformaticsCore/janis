from Pipeline.tool.commandtool import CommandTool, ToolInput
from Pipeline.types.common_data_types import Array, File


class Cat(CommandTool):

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