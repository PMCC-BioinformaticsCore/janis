from Pipeline.bioinformatics.tools.gatk4.gatk_4_0 import Gatk_4_0
from Pipeline.bioinformatics.tools.gatk4.markduplicates.base import Gatk4MarkDuplicatesBase


class Gatk4MarkDuplicates_4_0(Gatk_4_0, Gatk4MarkDuplicatesBase):
    pass


if __name__ == "__main__":
    print(Gatk4MarkDuplicates_4_0().help())
