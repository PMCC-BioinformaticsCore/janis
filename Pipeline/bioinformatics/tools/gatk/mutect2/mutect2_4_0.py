from Pipeline.bioinformatics.tools.gatk.gatk_4_0 import Gatk_4_0
from Pipeline.bioinformatics.tools.gatk.mutect2.base import GatkMutect2Base


class GatkMutect2_4_0(Gatk_4_0, GatkMutect2Base):
    pass


if __name__ == "__main__":
    print(GatkMutect2_4_0().help())