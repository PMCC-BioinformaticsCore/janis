from Pipeline.bioinformatics.tools.htslib.htslib_1_2_1 import HTSLib_1_2_1
from Pipeline.bioinformatics.tools.htslib.tabix.base import TabixBase


class Tabix_1_2_1(HTSLib_1_2_1, TabixBase):
    pass


if __name__ == "__main__":
    print(Tabix_1_2_1().help())
