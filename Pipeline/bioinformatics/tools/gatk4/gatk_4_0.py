from abc import ABC

class Gatk_4_0(ABC):

    @staticmethod
    def docker():
        return "broadinstitute/gatk:4.0.12.0"
