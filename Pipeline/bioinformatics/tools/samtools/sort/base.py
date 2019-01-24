from abc import ABC

from Pipeline import ToolInput, Filename, Int, String, Boolean, ToolOutput, Array
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.fasta import FastaWithDict
from Pipeline.bioinformatics.tools.samtools.samtoolstoolbase import SamToolsToolBase
from Pipeline.utils.metadata import ToolMetadata


class SamToolsSortBase(SamToolsToolBase, ABC):

    @staticmethod
    def tool():
        return "SamToolsSort"

    @classmethod
    def samtools_command(cls):
        return "sort"

    def inputs(self):
        return [
            *super(SamToolsSortBase, self).inputs(),
            *SamToolsSortBase.additional_inputs,

            ToolInput("bam", Bam(), position=10),

            ToolInput("outputFilename", Filename(extension=".bam"), position=5, prefix="-o",
                      doc="Output to FILE [stdout]."),
        ]

    def outputs(self):
        return [
            ToolOutput("output", Bam(), glob="$(inputs.outputFilename)"),
            ToolOutput("temporaryOutputs", Array(Bam(), optional=True), glob="*.tmp.*.bam",
                       doc="By default, any temporary files are written alongside the output file, "
                           "as out.bam.tmp.nnnn.bam, or if output is to standard output, "
                           "in the current directory as samtools.mmm.mmm.tmp.nnnn.bam.")
        ]

    def friendly_name(self):
        return "SamTools: Sort"

    def metadata(self):
        from datetime import date
        return ToolMetadata(
            creator="Michael Franklin",
            maintainer="Michael Franklin",
            maintainer_email="michael.franklin@petermac.org",
            date_created=date(2018, 12, 24),
            date_updated=date(2019, 1, 24),
            institution="Samtools",
            doi=None,
            citation=None, # find citation
            keywords=["samtools", "sort"],
            documentation_url="http://www.htslib.org/doc/samtools.html#DESCRIPTION",
            documentation="""Ensure SAMTOOLS.SORT is inheriting from parent metadata
    
---------------------------------------------------------------------------------------------------

Sort alignments by leftmost coordinates, or by read name when -n is used. An appropriate 
@HD-SO sort order header tag will be added or an existing one updated if necessary.

The sorted output is written to standard output by default, or to the specified file (out.bam) 
when -o is used. This command will also create temporary files tmpprefix.%d.bam as needed when 
the entire alignment data cannot fit into memory (as controlled via the -m option).

---------------------------------------------------------------------------------------------------

The following rules are used for ordering records.

If option -t is in use, records are first sorted by the value of the given alignment tag, and then 
by position or name (if using -n). For example, “-t RG” will make read group the primary sort key. 
The rules for ordering by tag are:

- Records that do not have the tag are sorted before ones that do.
- If the types of the tags are different, they will be sorted so that single character tags (type A) 
    come before array tags (type B), then string tags (types H and Z), then numeric tags (types f and i).
- Numeric tags (types f and i) are compared by value. Note that comparisons of floating-point values 
    are subject to issues of rounding and precision.
- String tags (types H and Z) are compared based on the binary contents of the tag using the C strcmp(3) function.
- Character tags (type A) are compared by binary character value.
- No attempt is made to compare tags of other types — notably type B array values will not be compared.

When the -n option is present, records are sorted by name. Names are compared so as to give a 
“natural” ordering — i.e. sections consisting of digits are compared numerically while all other 
sections are compared based on their binary representation. This means “a1” will come before 
“b1” and “a9” will come before “a10”. Records with the same name will be ordered according to 
the values of the READ1 and READ2 flags (see flags).

When the -n option is not present, reads are sorted by reference (according to the order of the 
@SQ header records), then by position in the reference, and then by the REVERSE flag.

*Note*

    Historically samtools sort also accepted a less flexible way of specifying the 
    final and temporary output filenames:
    
    |   samtools sort [-f] [-o] in.bam out.prefix
    
    This has now been removed. The previous out.prefix argument (and -f option, if any) 
    should be changed to an appropriate combination of -T PREFIX and -o FILE. The previous -o 
    option should be removed, as output defaults to standard output."""
        )

    @staticmethod
    def stdout():
        return "$(inputs.outputFilename)"

    additional_inputs = [
        ToolInput("compression", Int(optional=True), prefix="-l",
                  doc="Set the desired compression level for the final output file, ranging from "
                      "0 (uncompressed) or 1 (fastest but minimal compression) to 9 (best compression "
                      "but slowest to write), similarly to gzip(1)'s compression level setting.\n"
                      "If -l is not used, the default compression level will apply."),
        ToolInput("maximumMemory", String(optional=True), prefix="-m",
                  doc="Approximately the maximum required memory per thread, specified  either in bytes or "
                      "with a K, M, or G suffix [768 MiB]. "
                      "To prevent sort from creating a huge number of temporary files, "
                      "it enforces a minimum value of 1M for this setting."),
        ToolInput("sortByReadNames", Boolean(optional=True), prefix="-n",
                  doc="Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates."),
        ToolInput("outputType", String(optional=True), prefix="-O",
                  doc="Write the final output as sam, bam, or cram. By default, samtools tries "
                      "to select a format based on the -o filename extension; if output is to "
                      "standard output or no format can be deduced, bam is selected."),
        ToolInput("temporaryFilesPrefix", String(optional=True), prefix="-T",
                  doc="Write temporary files to PREFIX.nnnn.bam, or if the specified PREFIX is an "
                      "existing directory, to PREFIX/samtools.mmm.mmm.tmp.nnnn.bam, where mmm is unique "
                      "to this invocation of the sort command.\n"
                      "By default, any temporary files are written alongside the output file, "
                      "as out.bam.tmp.nnnn.bam, or if output is to standard output, in the current "
                      "directory as samtools.mmm.mmm.tmp.nnnn.bam."),
        ToolInput("threads", Int(optional=True), prefix="-@",
                  doc="Set number of sorting and compression threads. By default, operation is single-threaded.")
    ]
