from Pipeline.bioinformatics.tools.gatk4.gatk_4_0 import Gatk_4_0
from Pipeline.bioinformatics.tools.gatk4.mergesamfiles.base import Gatk4MergeSamFilesBase


class Gatk4MergeSamFiles_4_0(Gatk_4_0, Gatk4MergeSamFilesBase):
    pass


if __name__ == "__main__":
    print(Gatk4MergeSamFiles_4_0().help())
