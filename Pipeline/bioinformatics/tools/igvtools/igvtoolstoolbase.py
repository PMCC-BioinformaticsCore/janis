from abc import ABC, abstractmethod
from Pipeline import CommandTool


class IgvToolsToolBase(CommandTool, ABC):

    @classmethod
    @abstractmethod
    def igvtools_command(cls):
        raise Exception("Subclass must implement the igvtools_command method: expects one of: "
                        "[ toTDF, count, index, sort, version ]")

    @classmethod
    def base_command(cls):
        return ["samtools", cls.igvtools_command()]

    def doc(self):
        return """
    The igvtools utility provides a set of tools for pre-processing data files. 
    File names must contain an accepted file extension, e.g. test-xyz.bam. Tools include:

    toTDF: 
        Converts a sorted data input file to a binary tiled data (.tdf) file. Used to preprocess large datasets 
        for improved IGV performance. Supported input file formats: .cn, .gct, .igv, .res, .snp, and .wig
        Note: This tool was previously known as _tile_ 
    
    count:
        Computes average alignment or feature density for over a specified window size across the genome and 
        outputs a binary tiled data .tdf file, text .wig file, or both depending on inputs.
        Used to create a track that can be displayed in IGV, for example as a bar chart.
        Supported input file formats: .aligned, .bam, .bed, .psl, .pslx, and .sam
    
    index
        Creates an index file for an ASCII alignment or feature file.
        Index files are required for loading alignment files into IGV, and can significantly improve 
        performance for large feature files. The file must be in sorted order by start position.
        Supported input file formats: .aligned, .bed, .psl, .sam, .bam, and .vcf (v3.2)

    sort
        Sorts the input file by start position. 
        Used to prepare data files for tools that required sorted input files.
        Supported input file formats: .aligned, .bed, .cn, .igv, .psl, .sam, .bam, and .vcf

    Documentation: https://software.broadinstitute.org/software/igv/igvtools""".strip()

    @staticmethod
    @abstractmethod
    def docker():
        raise Exception("An error likely occurred when resolving the method order for docker for the igvtools classes "
                        "or you're trying to execute the docker method of the base class (ie, don't do that). "
                        "The method order resolution must preference Gatkbase subclasses, "
                        "and the subclass must contain a definition for docker.")

    def arguments(self):
        return []