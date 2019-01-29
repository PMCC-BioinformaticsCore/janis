from janis import CommandTool, String, ToolInput, ToolOutput, Stdout


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
            ToolInput("input", String())
        ]

    def outputs(self):
        return [ToolOutput("output", Stdout())]
