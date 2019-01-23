from abc import ABC
from typing import List

from Pipeline import CommandTool, ToolOutput, ToolInput, Boolean, Int, String, File
from Pipeline.bioinformatics.data_types.vcf import TabixIdx, VcfIdx, CompressedVcf


class TabixBase(CommandTool, ABC):
    @staticmethod
    def tool():
        return "tabix"

    @staticmethod
    def base_command():
        return "tabix"

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("file", CompressedVcf(), position=8,
                      doc="File from which to create the index. The input data file must be position sorted and "
                          "compressed by bgzip which has a gzip(1) like interface."),
            ToolInput("preset", String(optional=True), prefix="--preset", position=2,
                      doc="-p: Input format for indexing. Valid values are: gff, bed, sam, vcf. This option should "
                          "not be applied together with any of -s, -b, -e, -c and -0; it is not used for data "
                          "retrieval because this setting is stored in the index file. [gff]", default="vcf"),
            *self.additional_args
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("output", TabixIdx(), glob="$(inputs.file.basename)")
        ]

    @staticmethod
    def requirements():
        import cwlgen.cwlgen as cwl
        return [
            cwl.InitialWorkDirRequirement([
                cwl.InitialWorkDirRequirement.Dirent("$(inputs.file)")
            ])
        ]

    @staticmethod
    def docurl():
        return "http://www.htslib.org/doc/tabix.html"

    def doc(self):
        return """
    tabix – Generic indexer for TAB-delimited genome position files
    
    Tabix indexes a TAB-delimited genome position file in.tab.bgz and creates an index file (in.tab.bgz.tbi or 
    in.tab.bgz.csi) when region is absent from the command-line. The input data file must be position sorted 
    and compressed by bgzip which has a gzip(1) like interface.

    After indexing, tabix is able to quickly retrieve data lines overlapping regions specified in the format 
    "chr:beginPos-endPos". (Coordinates specified in this region format are 1-based and inclusive.)

    Fast data retrieval also works over network if URI is given as a file name and in this case the 
    index file will be downloaded if it is not present locally.""".strip()

    additional_args = [
        ToolInput("zeroBased", Boolean(optional=True), prefix="--zero-based", position=1,
                  doc="-0: Specify that the position in the data file is 0-based (e.g. UCSC files) rather than 1-based."),

        ToolInput("begin", Int(optional=True), prefix="--begin", position=4,
                  doc="-b: Column of start chromosomal position. [4]"),

        ToolInput("comment", String(optional=True), prefix="--comment", position=7,
                  doc="-c: Skip lines started with character CHAR. [#]"),

        ToolInput("csi", Boolean(optional=True), prefix="--csi", position=1,
                  doc="-C: Produce CSI format index instead of classical tabix or BAI style indices."),

        ToolInput("end", Int(optional=True), prefix="--end", position=5,
                  doc="-e: Column of end chromosomal position. The end column can be the same as the start column. [5]"),

        ToolInput("force", Boolean(optional=True), prefix="--force", position=1,
                  doc="-f: Force to overwrite the index file if it is present."),

        ToolInput("minShift", Int(optional=True), prefix="--min-shift", position=1,
                  doc="-m: set minimal interval size for CSI indices to 2^INT [14]"),

        ToolInput("sequence", Int(optional=True), prefix="--sequence", position=3,
                  doc="-s: Column of sequence name. Option -s, -b, -e, -S, -c and -0 are all stored "
                      "in the index file and thus not used in data retrieval. [1]"),

        ToolInput("skipLines", Int(optional=True), prefix="--skip-lines", position=6,
                  doc="-S: Skip first INT lines in the data file. [0]"),

        ToolInput("printHeader", Boolean(optional=True), prefix="--print-header", position=1,
                  doc="-h: Print also the header/meta lines."),

        ToolInput("onlyHeader", Boolean(optional=True), prefix="--only-header", position=1,
                  doc="-H: Print only the header/meta lines."),

        ToolInput("listChroms", Boolean(optional=True), prefix="--list-chroms", position=1,
                  doc="-l: List the sequence names stored in the index file."),

        ToolInput("reheader", File(optional=True), prefix="--reheader", position=1,
                  doc="-r: Replace the header with the content of FILE"),

        ToolInput("regions", File(optional=True), prefix="--regions", position=11,
                  doc="-R: Restrict to regions listed in the FILE. The FILE can be BED file "
                      "(requires .bed, .bed.gz, .bed.bgz file name extension) or a TAB-delimited "
                      "file with CHROM, POS, and, optionally, POS_TO columns, where positions are "
                      "1-based and inclusive. When this option is in use, the input file may not be sorted."),

        ToolInput("targets", File(optional=True), prefix="--targets", position=11,
                  doc="-T: Similar to -R but the entire input will be read "
                      "sequentially and regions not listed in FILE will be skipped"),
    ]
