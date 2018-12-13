from Pipeline import Tool, Input, String, ToolInput


class Echo(Tool):

    inp = ToolInput("inp", String())

    @staticmethod
    def tool():
        return "Echo"

    @staticmethod
    def base_command():
        return "echo"
