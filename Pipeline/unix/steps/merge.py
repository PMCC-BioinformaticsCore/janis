from Pipeline import File, Array, CommandTool, ToolInput, ToolOutput


class Merge(CommandTool):

    @staticmethod
    def tool():
        return "merge"

    @staticmethod
    def base_command():
        return ["cat"]

    def inputs(self):
        return [
            ToolInput("files", Array(File()))
        ]

    def outputs(self):
        return [
            ToolOutput("merged", File())
        ]
