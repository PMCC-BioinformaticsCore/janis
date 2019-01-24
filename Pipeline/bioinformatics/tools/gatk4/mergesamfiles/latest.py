from Pipeline.bioinformatics.tools.gatk4.gatk_latest import Gatk4Latest
from Pipeline.bioinformatics.tools.gatk4.mergesamfiles.base import Gatk4MergeSamFilesBase


class Gatk4MergeSamFilesLatest(Gatk4Latest, Gatk4MergeSamFilesBase):
    pass


if __name__ == "__main__":
    print(Gatk4MergeSamFilesLatest().help())
