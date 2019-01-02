from abc import ABC, abstractmethod

from Pipeline import CommandTool


class SamToolsToolBase(CommandTool, ABC):

    @staticmethod
    @classmethod
    def samtools_command(cls):
        raise Exception("Subclass must implement the samtools_command method: expects one of: ["
                        "   view, sort, index, idxstats, flagstat, stats, bedcov, depth, "
                        "   merge, faidx, fqidx, tview, split, quickcheck, dict, fixmate, "
                        "   mpileup, flags, fastq/a, collate, refheader, cat, rmdup, "
                        "   addreplacerg, calmd, targetcut, phase, depad, markdup"
                        "]")

    @classmethod
    def base_command(cls):
        return ["samtools", cls.samtools_command()]

    def inputs(self):
        return [

        ]

    @staticmethod
    def doc():
        return """
    Samtools is a set of utilities that manipulate alignments in the BAM format. It imports from 
    and exports to the SAM (Sequence Alignment/Map) format, does sorting, merging and indexing, 
    and allows to retrieve reads in any regions swiftly.

    Samtools is designed to work on a stream. It regards an input file `-' as the standard input (stdin) 
    and an output file `-' as the standard output (stdout). Several commands can thus be combined with 
    Unix pipes. Samtools always output warning and error messages to the standard error output (stderr).

    Samtools is also able to open a BAM (not SAM) file on a remote FTP or HTTP server if the BAM file 
    name starts with `ftp://' or `http://'. Samtools checks the current working directory for the index 
    file and will download the index upon absence. Samtools does not retrieve the entire alignment file 
    unless it is asked to do so.

    Documentation: http://www.htslib.org/doc/samtools.html#DESCRIPTION""".strip()

    @staticmethod
    @abstractmethod
    def docker():
        raise Exception("An error likely occurred when resolving the method order for docker for the samtools classes "
                        "or you're trying to execute the docker method of the base class (ie, don't do that). "
                        "The method order resolution must preference Gatkbase subclasses, "
                        "and the subclass must contain a definition for docker.")

    def arguments(self):
        return []