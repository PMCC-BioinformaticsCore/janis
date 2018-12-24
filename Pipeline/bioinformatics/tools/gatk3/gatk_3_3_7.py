from abc import ABC

from Pipeline.bioinformatics.tools.gatk3.gatk3base import Gatk3Base


class Gatk_3_3_7(Gatk3Base, ABC):

    @staticmethod
    def docker():
        return "broad/gatk:3."
