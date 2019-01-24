from Pipeline.bioinformatics.tools.gatk4.gatk_latest import Gatk4Latest
from Pipeline.bioinformatics.tools.gatk4.markduplicates.base import Gatk4MarkDuplicatesBase


class Gatk4MarkDuplicatesLatest(Gatk4Latest, Gatk4MarkDuplicatesBase):
    pass

if __name__ == "__main__":
    print(Gatk4MarkDuplicatesLatest().help())
