from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline import File, String, Int, Tool, ToolOutput, ToolInput


class PicardSortSam(Tool):
    inputFileName_sortSam = ToolInput("inputFileName_sortSam", File())
    outputFileName_sortSam = ToolInput("outputFileName_sortSam", String())
    tmpdir = ToolInput("tmpdir", String())
    validation_stringency = ToolInput("validation_stringency", String())
    maxRecordsInRam = ToolInput("maxRecordsInRam", Int(optional=True))

    out = ToolOutput("out", Bam()),             # Bam file
    indexes = ToolOutput("indexes", File())     # Bai Index

    @staticmethod
    def tool():
        return "picard-sortsam"

    @staticmethod
    def base_command():
        return "javac"
