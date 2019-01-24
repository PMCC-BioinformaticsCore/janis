from Pipeline import File, ToolInput
from Pipeline.bioinformatics.tools.gatk4.gatk_4_0 import Gatk_4_0
from Pipeline.bioinformatics.tools.gatk4.printreads.base import Gatk4PrintReadsBase


class Gatk4PrintReads_4_0(Gatk_4_0, Gatk4PrintReadsBase):
    pass

if __name__ == "__main__":
    print(Gatk4PrintReads_4_0().help())
