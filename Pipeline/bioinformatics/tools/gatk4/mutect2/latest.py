from Pipeline.bioinformatics.tools.gatk4.gatk_latest import Gatk4Latest
from Pipeline.bioinformatics.tools.gatk4.mutect2.base import Gatk4Mutect2Base


class Gatk4Mutect2Latest(Gatk4Latest, Gatk4Mutect2Base):
    pass


if __name__ == "__main__":
    print(Gatk4Mutect2Latest().help())
