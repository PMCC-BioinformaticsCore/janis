from Pipeline.bioinformatics.tools.gatk4.applybqsr.base import Gatk4ApplyBqsrBase
from Pipeline.bioinformatics.tools.gatk4.gatk_latest import Gatk4Latest


class Gatk4ApplyBqsrLatest(Gatk4Latest, Gatk4ApplyBqsrBase):
    pass


if __name__ == "__main__":
    print(Gatk4ApplyBqsrLatest().help())
