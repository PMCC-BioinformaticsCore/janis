from Pipeline import File, Array, CommandTool, ToolInput, ToolOutput


class Merge(CommandTool):
    files = ToolInput("files", Array(File()))
    merged = ToolOutput("merged", File())

    @staticmethod
    def tool():
        return "merge"

    @staticmethod
    def base_command():
        return ["cat"]