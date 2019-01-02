from Pipeline.bioinformatics.tools.gatk4.gatk_latest import Gatk4Latest
from Pipeline.bioinformatics.tools.gatk4.sortsam.base import Gatk4SortSamBase


class Gatk4SortSamLatest(Gatk4Latest, Gatk4SortSamBase):
    pass


if __name__ == "__main__":
    print(Gatk4SortSamLatest().help())
