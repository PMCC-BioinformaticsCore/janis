from Pipeline.bioinformatics.tools.samtools.samtoolslatest import SamToolsLatest
from Pipeline.bioinformatics.tools.samtools.view.base import SamToolsViewBase


class SamToolsViewLatest(SamToolsLatest, SamToolsViewBase):
    pass


if __name__ == "__main__":
    print(SamToolsViewLatest().help())
