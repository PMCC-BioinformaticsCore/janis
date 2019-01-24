from Pipeline.bioinformatics.tools.bwa.bwa_0_7_15 import Bwa_0_7_15
from Pipeline.bioinformatics.tools.bwa.mem.membase import BwaMemBase


class BwaMem_0_7_15(Bwa_0_7_15, BwaMemBase):
    pass


if __name__ == "__main__":
    print(BwaMem_0_7_15().help())
