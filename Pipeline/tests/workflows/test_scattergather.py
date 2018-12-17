from typing import List

from Pipeline import Workflow, Input, File, Directory, Int, Output, Array, String, Step, CommandTool, ToolInput, \
    ToolOutput
from Pipeline.tool.expressiontool import ExpressionTool


class Fasta(File):

    @staticmethod
    def name() -> str:
        return "Fasta"

    @staticmethod
    def doc() -> str:
        return "FASTA file"


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


class Merge(CommandTool):
    files = ToolInput("files", Array(File()))
    merged = ToolOutput("merged", File())

    @staticmethod
    def tool():
        return "merge"

    @staticmethod
    def base_command():
        return ["cat"]


class ParseFile(ExpressionTool):
    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("file", File())
        ]

    def outputs(self):
        return [
            ToolOutput("lines", Array(String()))
        ]


class Parse(ExpressionTool):
    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("line", String())
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("sampleName", String()),
            ToolOutput("fastq1", File()),
            ToolOutput("fastq2", File())
        ]


class TestScatterGather:


    def subworkflow(self):
        sw = Workflow("subworkflow")

        line = Input("line", String())
        blastDbNt = Input("blastDb", Directory())
        blastNumThreads = Input("blastNumThreads", Int(optional=True))

        contigd = Output("contigd", File())
        blastd = Output("blastd", File())

        step1 = Step("parse", Parse())
        step2 = Step("fq2fa", Fq2fa())
        step3 = Step("idba", Idba())
        step4 = Step("blast", Blast())

        sw.add_edges([
            (line, step1),

            (step1.sampleName, step2.sampleName),
            (step1.fastq1, step2.fastq1),
            (step1.fastq2, step2.fastq2),

            (step1.sampleName, step3.sampleName),
            (step2, step3),

            (step1.sampleName, step4),
            (step3, step4),
            (blastDbNt, step4),
            (blastNumThreads, step4),

            (step3, contigd),
            (step4, blastd)
        ])

        return sw


    def workflow(self):
        w = Workflow("Parkvile-data-workflow-demo")
        inputSamplesFile = Input("inputSamplesFile", File())
        blastDbNt = Input("blastDb", Directory())
        blastNumThreads = Input("blastNumThreads", Int(optional=True))

        contigd = Output("contigd", Array(File()))
        blastd = Output("blastd", Array(File()))
        merged = Output("merged", File())

        step1 = Step("parse-lines", ParseFile())
        step2 = Step("subworkflow", self.subworkflow())
        step3 = Step("merge", Merge())

        w.add_edges([
            (inputSamplesFile, step1),

            (step1.lines, step2.line),
            (blastDbNt, step2),
            (blastNumThreads, step2),

            (step2.blastd, step3.files),

            (step2.blastd, blastd),
            (step2.contigd, contigd),
            (step3.files, merged)
        ])

        return w



if __name__ == "__main__":
    w = TestScatterGather().workflow()
    w.dump_cwl(to_disk=True)
    print(w.cwl())