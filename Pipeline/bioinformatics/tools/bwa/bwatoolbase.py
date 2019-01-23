from abc import ABC, abstractmethod

from Pipeline import CommandTool, ToolInput, ToolArgument


class BwaToolBase(CommandTool, ABC):

    @staticmethod
    def base_command():
        return ["bwa"]

    def inputs(self):
        return [

        ]

    def doc(self):
        return """
    BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human 
    genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for 
    Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. 
    BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is 
    the latest, is generally recommended for high-quality queries as it is faster and more accurate. 
    BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.
    
    Documentation: http://bio-bwa.sourceforge.net/bwa.shtml
        """.strip()

    @staticmethod
    @abstractmethod
    def docker():
        raise Exception("An error likely occurred when resolving the method order for docker for the Gatk classes "
                        "or you're trying to execute the docker method of the base class (ie, don't do that). "
                        "The method order resolution must preference Gatkbase subclasses, "
                        "and the subclass must contain a definition for docker.")

    def arguments(self):
        return []
