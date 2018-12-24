from Pipeline.bioinformatics.tools.gatk.gatk_4_0 import Gatk_4_0
from Pipeline.bioinformatics.tools.gatk.baserecalibrator.base import GatkRecalibratorBase


class GatkBaseRecalibrator_4_0(Gatk_4_0, GatkRecalibratorBase):
    pass


if __name__ == "__main__":
    print(GatkBaseRecalibrator_4_0().help())
