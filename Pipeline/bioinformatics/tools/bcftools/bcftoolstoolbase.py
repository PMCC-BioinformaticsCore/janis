from abc import ABC, abstractmethod
from Pipeline import CommandTool


class BcfToolsToolBase(CommandTool, ABC):

    @classmethod
    @abstractmethod
    def bcftools_command(cls):
        raise Exception("Subclass must implement the bcftools_command method: expects one of: ["
                        "   annotate, call, cnv, concat, consensus, convert, csq, "
                        "   filter, gtcheck, index, isec, merge, mpileup, norm, "
                        "   plugin, polysomy, query, reheader, roh, sort, stats, view"
                        "]")

    @classmethod
    def base_command(cls):
        return ["bcftools", cls.bcftools_command()]

    def doc(self):
        return """
    BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary 
    counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed.

    Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming 
    from a pipe. \Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will 
    work in most, but not all situations. In general, whenever multiple VCFs are read simultaneously, 
    they must be indexed and therefore also compressed.
    
    BCFtools is designed to work on a stream. It regards an input file "-" as the standard input (stdin) 
    and outputs to the standard output (stdout). Several commands can thus be combined with Unix pipes.

    Documentation: https://samtools.github.io/bcftools/bcftools.html""".strip()

    @staticmethod
    @abstractmethod
    def docker():
        raise Exception("An error likely occurred when resolving the method order for docker for the bcftools classes "
                        "or you're trying to execute the docker method of the base class (ie, don't do that). "
                        "The method order resolution must preference Gatkbase subclasses, "
                        "and the subclass must contain a definition for docker.")

    def arguments(self):
        return []
