from Pipeline.bioinformatics.tools.samtools.samtoolstoolbase import SamToolsToolBase


class SamTools_1_7(SamToolsToolBase):
    @staticmethod
    def docker():
        return "biocontainers/samtools:v1.7.0_cv3"
