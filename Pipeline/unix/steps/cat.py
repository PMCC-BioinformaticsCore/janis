from Pipeline.tool.tool import Tool, ToolInput
from Pipeline.types.common_data_types import Array, File


class Cat(Tool):

    files = ToolInput("files", Array(File()))

    @staticmethod
    def tool():
        return "CAT"

    @staticmethod
    def base_command():
        return "cat"
