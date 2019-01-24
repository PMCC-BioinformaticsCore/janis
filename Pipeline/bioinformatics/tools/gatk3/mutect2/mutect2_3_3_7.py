from Pipeline.bioinformatics.tools.gatk3.gatk_3_3_7 import Gatk_3_3_7
from Pipeline.bioinformatics.tools.gatk3.mutect2.base import Gatk3Mutect2Base


class Gatk3Mutect2_3_3_7(Gatk_3_3_7, Gatk3Mutect2Base):
    pass


if __name__ == "__main__":
    print(Gatk3Mutect2_3_3_7().help())
