from Pipeline.bioinformatics.tools.gatk4.gatk_4_0 import Gatk_4_0
from Pipeline.bioinformatics.tools.gatk4.mutect2.base import Gatk4Mutect2Base


class GatkMutect2_4_0(Gatk_4_0, Gatk4Mutect2Base):
    pass


if __name__ == "__main__":
    print(GatkMutect2_4_0().help())
