from Pipeline.bioinformatics.tools.bwa.bwabase import BwaBase


class Bwa_0_7_15(BwaBase):
    @staticmethod
    def docker():
        return "biocontainers/bwa:v0.7.15_cv3"
