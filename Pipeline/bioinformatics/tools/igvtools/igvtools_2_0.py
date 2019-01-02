from Pipeline.bioinformatics.tools.igvtools.igvtoolsbase import IgvToolsBase


class IgvTools_2_0(IgvToolsBase):
    @staticmethod
    def docker():
        return "maxulysse/igvtools:2.0.0"
