from Pipeline import CommandTool, ToolInput, File, String, ToolOutput


class Idba(CommandTool):

    sampleName = ToolInput("sampleName", String())
    idbaInputFasta = ToolInput("idbaInputFasta", File())

    idbaOutputFasta = ToolOutput("idbaOutputFasta", File())

    @staticmethod
    def tool():
        return "idba"

    @staticmethod
    def base_command():
        return "blast"