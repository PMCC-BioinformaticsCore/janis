from Pipeline.bioinformatics.tools.gatk4.applybqsr.base import Gatk4ApplyBqsrBase
from Pipeline.bioinformatics.tools.gatk4.gatk_4_0 import Gatk_4_0


class Gatk4ApplyBqsr_4_0(Gatk_4_0, Gatk4ApplyBqsrBase):
    pass


if __name__ == "__main__":
    print(Gatk4ApplyBqsr_4_0().help())
