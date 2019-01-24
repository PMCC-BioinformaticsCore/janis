from Pipeline.bioinformatics.tools.gatk4.gatk_4_0 import Gatk_4_0
from Pipeline.bioinformatics.tools.gatk4.genotypeconcordance.base import Gatk4GenotypeConcordanceBase


class Gatk4GenotypeConcordance_4_0(Gatk_4_0, Gatk4GenotypeConcordanceBase):
    pass


if __name__ == "__main__":
    print(Gatk4GenotypeConcordance_4_0().help())
