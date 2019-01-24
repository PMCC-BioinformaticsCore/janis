from Pipeline.bioinformatics.tools.gatk3.gatk3_latest import Gatk3Latest
from Pipeline.bioinformatics.tools.gatk3.mutect2.base import Gatk3Mutect2Base


class Gatk3Mutect2Latest(Gatk3Latest, Gatk3Mutect2Base):
    pass


if __name__ == "__main__":
    print(Gatk3Mutect2Latest().help())
