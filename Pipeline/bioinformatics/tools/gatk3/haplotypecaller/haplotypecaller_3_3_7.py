from Pipeline.bioinformatics.tools.gatk.gatk_4_0 import Gatk_4_0
from Pipeline.bioinformatics.tools.gatk.haplotypecaller.base import GatkHaplotypeCallerBase


class GatkHaplotypeCaller_4_0(Gatk_4_0, GatkHaplotypeCallerBase):
    pass


if __name__ == "__main__":
    print(GatkHaplotypeCaller_4_0().help())
