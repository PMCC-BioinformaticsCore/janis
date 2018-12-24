from Pipeline.bioinformatics.tools.gatk.gatk_latest import GatkLatest
from Pipeline.bioinformatics.tools.gatk.mutect2.base import GatkMutect2Base


class GatkMutect2Latest(GatkLatest, GatkMutect2Base):
    pass


if __name__ == "__main__":
    print(GatkMutect2Latest().help())
