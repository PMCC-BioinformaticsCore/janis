from Pipeline.bioinformatics.tools.gatk4.gatk_latest import Gatk4Latest
from Pipeline.bioinformatics.tools.gatk4.printreads.base import Gatk4PrintReadsBase


class Gatk4PrintReadsLatest(Gatk4Latest, Gatk4PrintReadsBase):
    pass


if __name__ == "__main__":
    print(Gatk4PrintReadsLatest().help())
