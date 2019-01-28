from Pipeline import CommandTool, Input, String, ToolInput


class Echo(CommandTool):

    @staticmethod
    def tool():
        return "Echo"

    def friendly_name(self):
        return "Echo - Print to console"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [
            ToolInput("inp", String())
        ]

    def outputs(self):
        return []
