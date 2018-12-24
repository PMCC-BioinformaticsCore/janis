from Pipeline.bioinformatics.tools.gatk4.gatk_4_0 import Gatk_4_0
from Pipeline.bioinformatics.tools.gatk4.baserecalibrator.base import Gatk4BaseRecalibratorBase


class Gatk4BaseRecalibrator_4_0(Gatk_4_0, Gatk4BaseRecalibratorBase):
    pass


if __name__ == "__main__":
    print(Gatk4BaseRecalibrator_4_0().help())
