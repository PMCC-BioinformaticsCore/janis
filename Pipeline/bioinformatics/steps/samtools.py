from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.sam import Sam
from Pipeline import String, CommandTool, ToolOutput, ToolInput, ToolArgument, Boolean, Int, File, Float, Array
from Pipeline.types.filename import Filename


class SamTools(CommandTool):
    inp = ToolInput("input", Bam(), doc="Input bam file.")
    out = ToolOutput("out", Sam(), glob="$(inputs.output_name")

    # Optional params
    outputName = ToolInput("outputName", Filename(extension="sam"), position=2, prefix="-o")

    @staticmethod
    def tool():
        return "samtools-view"

    @staticmethod
    def base_command():
        return ["samtools", "view"]

    @staticmethod
    def docker():
        return "biocontainers/samtools"

    def arguments(self):
        return [
            ToolArgument("-S", position=2),
            ToolArgument("-h", position=3),
            ToolArgument("-b", position=4)
        ]

    isBam = ToolInput("isBam", Boolean(), default=False, position=2, prefix="-b",
                      doc="output in BAM format")

    readsWithoutBits = ToolInput("readsWithoutBits", Int(optional=True), position=1, prefix="-F",
                                 doc="only include reads with none of the bits set in INT set in FLAG [0]")

    collapsecigar = ToolInput("collapsecigar", Boolean(), default=False, position=1, prefix="-B",
                              doc="collapse the backward CIGAR operation")

    readsingroup = ToolInput("readsingroup", String(optional=True), position=1,
                             prefix=" -r",
                             doc="only include reads in read group STR [null]")

    bedoverlap = ToolInput("bedoverlap", File(optional=True), position=1,
                           prefix=" -L",
                           doc="only include reads overlapping this BED FILE [null]")

    uncompressed = ToolInput("uncompressed", Boolean(), default=False, position=1, prefix=" -u",
                             doc="uncompressed BAM output (implies -b)")

    readtagtostrip = ToolInput("readtagtostrip", Array(String(), optional=True), position=1,
                               doc="read tag to strip (repeatable) [null]")

    readsquality = ToolInput("readsquality", Int(optional=True), position=1, prefix=" -q",
                             doc="only include reads with mapping quality >= INT [0]")

    readswithbits = ToolInput("readswithbits", Int(optional=True), position=1, prefix=" -f",
                              doc="only include reads with all bits set in INT set in FLAG [0]")

    cigar = ToolInput("cigar", Int(optional=True), position=1, prefix=" -m",
                      doc="only include reads with number of CIGAR operations consuming query sequence >= INT [0]")

    iscram = ToolInput("iscram", Boolean(), default=False, position=2, prefix=" -C",
                       doc="output in CRAM format")

    threads = ToolInput("threads", Int(optional=True), position=1, prefix=" -@",
                        doc="number of BAM compression threads [0]")

    fastcompression = ToolInput("fastcompression", Boolean(), default=False, position=1, prefix="-1",
                                doc="use fast BAM compression (implies -b)")

    samheader = ToolInput("samheader", Boolean(), default=False, position=1, prefix=" -h",
                          doc="include header in SAM output")

    count = ToolInput("count", Boolean(), default=False, position=1, prefix=" -c",
                      doc="prInt only the count of matching records")

    randomseed = ToolInput("randomseed", Float(optional=True), position=1, prefix=" -s",
                           doc="Integer part sets seed of random number generator [0]; sets fraction of templates to subsample [no subsampling]")

    referencefasta = ToolInput("referencefasta", File(optional=True), position=1, prefix=" -T",
                               doc="reference sequence FASTA FILE [null]")

    region = ToolInput("region", String(optional=True), position=5, doc="[region ...]")

    readsingroupfile = ToolInput("readsingroupfile", File(optional=True), position=1, prefix=" -R",
                                 doc="only include reads with read group listed in FILE [null]")

    readsinlibrary = ToolInput("readsinlibrary", String(optional=True), position=1, prefix=" -l",
                               doc="only include reads in library STR [null]")


if __name__ == "__main__":
    print(SamTools().help())
