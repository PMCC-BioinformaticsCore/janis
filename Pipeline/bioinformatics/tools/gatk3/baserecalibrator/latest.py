from Pipeline.bioinformatics.tools.gatk.gatk_latest import GatkLatest
from Pipeline.bioinformatics.tools.gatk.baserecalibrator.base import GatkRecalibratorBase


class GatkRecalibratorLatest(GatkLatest, GatkRecalibratorBase):
    pass

if __name__ == "__main__":
    print(GatkRecalibratorLatest().help())
