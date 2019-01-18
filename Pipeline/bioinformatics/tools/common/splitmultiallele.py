from typing import List

from Pipeline import CommandTool, ToolOutput, ToolInput, Filename
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.vcf import Vcf


class SplitMultiAllele(CommandTool):
    @staticmethod
    def tool():
        return "SplitMultiAllele"

    @staticmethod
    def base_command():
        return "vcfsplitmultiallele.sh"

    @staticmethod
    def docker():
        return "SEE (mfranklin's notes in) DOCUMENTATION"

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("input", Vcf()),
            ToolInput("outputFilename", Filename(extension=".vcf")),
            # ToolInput("reference", FastaWithDict()) # eventually
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("output", Vcf())
        ]

    @staticmethod
    def doc():
        return """
    VcfSplitMultiAllele.sh
    
    Currently stored at: '/researchers/jiaan.yu/WGS_pipeline/VcfSplitMultiAllele.sh’
    
    MFRANKLIN's notes:
        - At the moment, a random shell script is not portable as it exists in a location. I'd be happy
            with putting it in a docker container for when we work out how to run them on the cluster.
            
        - Some of the references are hard-coded, ie: the human reference genome. This won't correctly connect
            when using the execution engine, as I can't reference secondary (index) files and hence the execution
            engine won't correctly localize these.    
        """.strip()


if __name__ == "__main__":
    print(SplitMultiAllele().help())
