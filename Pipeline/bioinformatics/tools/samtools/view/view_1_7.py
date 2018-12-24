from Pipeline.bioinformatics.tools.samtools.samtools_1_7 import SamTools_1_7
from Pipeline.bioinformatics.tools.samtools.view.viewbase import SamToolsViewBase


class SamToolsView_1_7(SamTools_1_7, SamToolsViewBase):
    pass


if __name__ == "__main__":
    print(SamToolsView_1_7().help())
