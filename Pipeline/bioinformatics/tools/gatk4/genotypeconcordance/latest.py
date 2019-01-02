from Pipeline.bioinformatics.tools.gatk4.gatk_latest import Gatk4Latest
from Pipeline.bioinformatics.tools.gatk4.genotypeconcordance.base import Gatk4GenotypeConcordanceBase


class Gatk4GenotypeConcordanceLatest(Gatk4Latest, Gatk4GenotypeConcordanceBase):
    pass


if __name__ == "__main__":
    print(Gatk4GenotypeConcordanceLatest().help())
