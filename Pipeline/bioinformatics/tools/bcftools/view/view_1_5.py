from Pipeline.bioinformatics.tools.bcftools.bcftools_1_5 import BcfTools_1_5
from Pipeline.bioinformatics.tools.bcftools.view.base import BcfToolsViewBase


class BcfToolsView_1_5(BcfTools_1_5, BcfToolsViewBase):
    pass


if __name__ == "__main__":
    print(BcfToolsView_1_5().help())
