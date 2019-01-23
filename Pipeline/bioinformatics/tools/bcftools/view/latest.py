from Pipeline.bioinformatics.tools.bcftools.latest import BcfToolsLatest
from Pipeline.bioinformatics.tools.bcftools.view.base import BcfToolsViewBase


class BcfToolsViewLatest(BcfToolsLatest, BcfToolsViewBase):
    pass


if __name__ == "__main__":
    print(BcfToolsViewLatest().help())
