from abc import ABC

from Pipeline import ToolInput, ToolOutput, ToolArgument, Boolean, File, String, Int, Float, Filename
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.sam import Sam
from Pipeline.bioinformatics.tools.samtools.samtoolstoolbase import SamToolsToolBase

class SamToolsViewBase(SamToolsToolBase, ABC):
    @staticmethod
    def tool():
        return "SamToolsView"

    @staticmethod
    def base_command():
        # lol python: https://stackoverflow.com/a/26807879
        bc = super(SamToolsViewBase, SamToolsViewBase).base_command()
        if isinstance(bc, str): bc = [bc]
        return [*bc, "view"]

    def inputs(self):
        return [
            *super(SamToolsViewBase, self).inputs(),
            *SamToolsViewBase.additional_inputs,

            ToolInput("sam", Bam(), position=10),

            ToolInput("reference", FastaWithDict(optional=True), position=5, prefix="-T",
                      doc="A FASTA format reference FILE, optionally compressed by bgzip and ideally indexed "
                          "by samtools faidx. If an index is not present, one will be generated for you."),

            ToolInput("output", Filename(extension=".bam"), position=5, prefix="-o", doc="Output to FILE [stdout]."),

        ]

    @staticmethod
    def stdout():
        return "$(inputs.outputFilename)"

    def outputs(self):
        return [
            ToolOutput("out", Bam(), glob="$(inputs.outputFilename)"),
        ]

    @staticmethod
    def doc():
        return """
    SAMTOOLS: view
    
    With no options or regions specified, prints all alignments in the specified input alignment file 
    (in SAM, BAM, or CRAM format) to standard output in SAM format (with no header).
    
    You may specify one or more space-separated region specifications after the input filename to 
    restrict output to only those alignments which overlap the specified region(s). 
    Use of region specifications requires a coordinate-sorted and indexed input file (in BAM or CRAM format).
       
       
    Documentation: http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS
           """.strip()

    def arguments(self):
        return [
            ToolArgument("-S", position=2,
                         doc="Ignored for compatibility with previous samtools versions. Previously this option was "
                             "required if input was in SAM format, but now the correct format is automatically "
                             "detected by examining the first few characters of input."),
            ToolArgument("-h", position=3, doc="Include the header in the output."),
            ToolArgument("-b", position=4, doc="Output in the BAM format.")

        ]

    additional_inputs = [

        ToolInput("cramOutput", Boolean(optional=True), position=5, prefix="-C",
                  doc="Output in the CRAM format (requires -T)."),
        ToolInput("compressedBam", Boolean(optional=True), position=5, prefix="-1",
                  doc="Enable fast BAM compression (implies -b)."),
        ToolInput("uncompressedBam", Boolean(optional=True), position=5, prefix="-u",
                  doc="Output uncompressed BAM. This option saves time spent on compression/decompression and is "
                      "thus preferred when the output is piped to another samtools command."),
        ToolInput("onlyOutputHeader", Boolean(optional=True), position=5, prefix="-H", doc="Output the header only."),
        ToolInput("countAlignments", Boolean(optional=True), position=5, prefix="-c",
                  doc="Instead of printing the alignments, only count them and print the total number. "
                      "All filter options, such as -f, -F, and -q, are taken into account."),
        # ToolInput("", Boolean(), position=5, prefix="-?", doc="Output long help and exit immediately."),
        ToolInput("writeAlignments", File(optional=True), position=5, prefix="-U",
                  doc="Write alignments that are not selected by the various filter options to FILE. "
                      "When this option is used, all alignments (or all alignments intersecting the regions specified) "
                      "are written to either the output file or this file, but never both."),
        ToolInput("inputTSV", File(optional=True), position=5, prefix="-t",
                  doc="A tab-delimited FILE. Each line must contain the reference name in the first column and the "
                      "length of the reference in the second column, with one line for each distinct reference. "
                      "Any additional fields beyond the second column are ignored. This file also defines the order "
                      "of the reference sequences in sorting. If you run: `samtools faidx <ref.fa>', the resulting "
                      "index file <ref.fa>.fai can be used as this FILE."),
        ToolInput("onlyOverlapping", File(optional=True), position=5, prefix="-L",
                  doc="Only output alignments overlapping the input BED FILE [null]."),
        ToolInput("useMultiRegionIterator", Boolean(optional=True), position=5, prefix="-M",
                  doc="Use the multi-region iterator on the union of the BED file and command-line region arguments. "
                      "This avoids re-reading the same regions of files so can sometimes be much faster. "
                      "Note this also removes duplicate sequences. Without this a sequence that overlaps multiple "
                      "regions specified on the command line will be reported multiple times."),
        ToolInput("outputAlignmentsInReadGroup", String(optional=True), position=5, prefix="-r",
                  doc="Output alignments in read group STR [null]. Note that records with no RG tag will also be "
                      "output when using this option. This behaviour may change in a future release."),
        ToolInput("outputAlignmentsInFileReadGroups", File(optional=True), position=5, prefix="-R",
                  doc="Output alignments in read groups listed in FILE [null]. Note that records with no RG tag "
                      "will also be output when using this option. This behaviour may change in a future release."),
        ToolInput("mapqThreshold", Int(optional=True), position=5, prefix="-q", doc="Skip alignments with MAPQ smaller than INT [0]."),
        ToolInput("outputAlignmentsInLibrary", String(optional=True), position=5, prefix="-l", doc="Only output alignments in library STR [null]."),
        ToolInput("outputAlignmentsMeetingCIGARThreshold", Int(optional=True), position=5, prefix="-m",
                  doc="Only output alignments with number of CIGAR bases consuming query sequence ≥ INT [0]"),
        ToolInput("outputAlignmentsWithBitsSet", Int(optional=True), position=5, prefix="-f",
                  doc="Only output alignments with all bits set in INT present in the FLAG field. "
                      "INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or "
                      "in octal by beginning with `0' (i.e. /^0[0-7]+/) [0]."),
        ToolInput("doNotOutputAlignmentsWithBitsSet", Int(optional=True), position=5, prefix="-F",
                  doc="Do not output alignments with any bits set in INT present in the FLAG field. "
                      "INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or "
                      "in octal by beginning with `0' (i.e. /^0[0-7]+/) [0]."),
        ToolInput("doNotOutputAlignmentsWithAllBitsSet", Int(optional=True), position=5, prefix="-G",
                  doc="Do not output alignments with all bits set in INT present in the FLAG field. "
                      "This is the opposite of -f such that -f12 -G12 is the same as no filtering at all. "
                      "INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or "
                      "in octal by beginning with `0' (i.e. /^0[0-7]+/) [0]."),
        ToolInput("readTagToExclude", String(optional=True), position=5, prefix="-x", doc="Read tag to exclude from output (repeatable) [null]"),
        ToolInput("collapseBackwardCIGAR", Boolean(optional=True), position=5, prefix="-B", doc="Collapse the backward CIGAR operation."),
        ToolInput("subsamplingProportion", Float(optional=True), position=5, prefix="-s",
                  doc="Output only a proportion of the input alignments. This subsampling acts in the same "
                      "way on all of the alignment records in the same template or read pair, so it never "
                      "keeps a read but not its mate. The integer and fractional parts of the -s INT.FRAC "
                      "option are used separately: the part after the decimal point sets the fraction of "
                      "templates/pairs to be kept, while the integer part is used as a seed that influences "
                      "which subset of reads is kept."),
        ToolInput("threads", Int(optional=True), position=5, prefix="-@",
                  doc="Number of BAM compression threads to use in addition to main thread [0].")

    ]
