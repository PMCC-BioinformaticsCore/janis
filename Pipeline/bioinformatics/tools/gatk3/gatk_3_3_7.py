from abc import ABC

# from Pipeline.bioinformatics.tools.gatk3.gatk3toolbase import Gatk3ToolBaseBase


class Gatk_3_3_7(ABC):

    @staticmethod
    def docker():
        return "broad/gatk:3."
