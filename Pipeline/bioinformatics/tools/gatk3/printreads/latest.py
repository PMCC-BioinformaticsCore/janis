from Pipeline.bioinformatics.tools.gatk3.gatk3_latest import Gatk3Latest
from Pipeline.bioinformatics.tools.gatk3.printreads.base import Gatk3PrintReadsBase


class Gatk3PrintReadsLatest(Gatk3Latest, Gatk3PrintReadsBase):
    pass


if __name__ == "__main__":
    print(Gatk3PrintReadsLatest().help())
