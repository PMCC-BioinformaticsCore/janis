from Pipeline import CommandTool, Input, String, ToolInput


class Echo(CommandTool):

    inp = ToolInput("inp", String())

    @staticmethod
    def tool():
        return "Echo"

    @staticmethod
    def base_command():
        return "echo"
