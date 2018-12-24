from Pipeline.bioinformatics.tools.gatk.gatk_latest import GatkLatest
from Pipeline.bioinformatics.tools.gatk.printreads.base import GatkPrintReadsBase


class GatkPrintReadsLatest(GatkLatest, GatkPrintReadsBase):
    pass

if __name__ == "__main__":
    print(GatkPrintReadsLatest().help())