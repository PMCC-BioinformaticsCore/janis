from Pipeline.bioinformatics.tools.gatk.gatk_latest import GatkLatest
from Pipeline.bioinformatics.tools.gatk.haplotypecaller.base import GatkHaplotypeCallerBase


class GatkHaplotypeCallerLatest(GatkLatest, GatkHaplotypeCallerBase):
    pass


if __name__ == "__main__":
    print(GatkHaplotypeCallerLatest().help())
