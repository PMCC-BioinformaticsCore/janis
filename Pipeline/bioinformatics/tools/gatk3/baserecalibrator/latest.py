from Pipeline.bioinformatics.tools.gatk3.gatk3_latest import Gatk3Latest
from Pipeline.bioinformatics.tools.gatk3.baserecalibrator.base import Gatk3RecalibratorBase


class Gatk3RecalibratorLatest(Gatk3Latest, Gatk3RecalibratorBase):
    pass


if __name__ == "__main__":
    print(Gatk3RecalibratorLatest().help())
