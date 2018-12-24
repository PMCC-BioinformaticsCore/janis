from Pipeline.bioinformatics.tools.samtools.samtoolsbase import SamToolsBase


class SamTools_1_7(SamToolsBase):
    @staticmethod
    def docker():
        return "biocontainers/samtools:v1.7.0_cv3"
