from Pipeline import CommandTool, ToolInput, File, String, ToolOutput, Directory, Int


class Blast(CommandTool):
    sampleName = ToolInput("sampleName", String())
    db = ToolInput("db", Directory())
    assembledContigs = ToolInput("assembledContigs", File())
    numThreads = ToolInput("numThreads", Int(optional=True))

    blastOutput = ToolOutput("blastOutput", File())

    @staticmethod
    def tool():
        return "blast"

    @staticmethod
    def base_command():
        return ["blast"]