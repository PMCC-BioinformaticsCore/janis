from Pipeline import CommandTool, Input, String, ToolInput


class Echo(CommandTool):

    @staticmethod
    def tool():
        return "Echo"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [
            ToolInput("inp", String())
        ]

    def outputs(self):
        return []
