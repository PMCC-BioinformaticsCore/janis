from abc import ABC

from Pipeline.bioinformatics.tools.gatk.gatk4base import Gatk4Base


class Gatk_4_0(Gatk4Base, ABC):

    @staticmethod
    def docker():
        return "broadinstitute/gatk:4.0.12.0"
