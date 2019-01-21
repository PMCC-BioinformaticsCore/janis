from Pipeline.bioinformatics.tools.bcftools.annotate.base import BcfToolsAnnotateBase
from Pipeline.bioinformatics.tools.bcftools.bcftools_1_5 import BcfTools_1_5


class BcfToolsAnnotate_1_5(BcfTools_1_5, BcfToolsAnnotateBase):
    pass


if __name__ == "__main__":
    print(BcfToolsAnnotate_1_5().help())
