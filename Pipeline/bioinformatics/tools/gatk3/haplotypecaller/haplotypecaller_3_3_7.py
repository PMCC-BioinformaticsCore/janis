from Pipeline.bioinformatics.tools.gatk3.gatk_3_3_7 import Gatk_3_3_7
from Pipeline.bioinformatics.tools.gatk3.haplotypecaller.base import Gatk3HaplotypeCallerBase


class GatkHaplotypeCaller_3_3_7(Gatk_3_3_7, Gatk3HaplotypeCallerBase):
    pass


if __name__ == "__main__":
    print(GatkHaplotypeCaller_3_3_7().help())
