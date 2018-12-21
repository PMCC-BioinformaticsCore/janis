from Pipeline import CommandTool, ToolInput, File, String, ToolOutput


class Fq2fa(CommandTool):

    sampleName = ToolInput("sampleName", String())
    fastq1 = ToolInput("fastq1", File())
    fastq2 = ToolInput("fastq2", File())

    outputFasta = ToolOutput("outputFasta", File())

    @staticmethod
    def tool():
        return "fq2fa"

    @staticmethod
    def base_command():
        return ["fq2fa"]