from typing import List

from janis import Workflow, Input, File, Directory, Int, Output, Array, String, Step, ToolInput, \
    ToolOutput
from bioinformatics import Blast
from bioinformatics import Fq2fa
from bioinformatics import Idba
from janis.tool.expressiontool import ExpressionTool
from janis.unix.tools.merge import Merge
from janis.unix.tools.scatterlines import ScatterLines


class ParseSamplesFile(ExpressionTool):

    @staticmethod
    def tool():
        return "parse"

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

    def expression(self):
        return """${
var l = inputs.line.split("\t");
return {
    sampleName: l[0], 
    fastq1: { class: "File", location: "file://" + l[1] }, 
    fastq2: { class: "File", location: "file://" + l[2] }
};
        }""".strip()


class TestScatterGather:

    def subworkflow(self):
        sw = Workflow("subworkflow")

        line = Input("line", String())
        blastDbNt = Input("blastDb", Directory())
        blastNumThreads = Input("blastNumThreads", Int(optional=True))

        contigd = Output("contigd", File())
        blastd = Output("blastd", File())

        step1 = Step("parse", ScatterLines())
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

        step1 = Step("parse-lines", ParseSamplesFile())
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