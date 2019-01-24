from Pipeline.bioinformatics.tools.gatk4.gatk_latest import Gatk4Latest
from Pipeline.bioinformatics.tools.gatk4.haplotypecaller.base import Gatk4HaplotypeCallerBase


class Gatk4HaplotypeCallerLatest(Gatk4Latest, Gatk4HaplotypeCallerBase):
    pass


if __name__ == "__main__":
    print(Gatk4HaplotypeCallerLatest().help())
