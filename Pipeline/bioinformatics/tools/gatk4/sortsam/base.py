from Pipeline import ToolInput, Filename, String, ToolArgument, Array, File, Int, Boolean, ToolOutput, Directory
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.fasta import FastaWithDict
from Pipeline.bioinformatics.data_types.sam import Sam
from Pipeline.bioinformatics.tools.gatk4.gatk4toolbase import Gatk4ToolBase


class Gatk4SortSamBase(Gatk4ToolBase):
    @classmethod
    def gatk_command(cls):
        return "SortSam"

    @staticmethod
    def tool():
        return "gatksortsam"

    @staticmethod
    def doc():
        return """
    Sorts a SAM/BAM/CRAM file.    
    
    Documentation: https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.3/org_broadinstitute_hellbender_tools_picard_sam_SortSam.php
    """.strip()

    def inputs(self):
        return [
            *super(Gatk4SortSamBase, self).inputs(),
            ToolInput("input", Bam(), prefix="-I", doc="The SAM/BAM/CRAM file to sort.", position=10),
            ToolInput("outputFilename", Filename(extension=".bam"), position=10, prefix="-O",
                      doc="The sorted SAM/BAM/CRAM output file."),
            ToolInput("sortOrder", String(), prefix="-SO", position=10,
                      doc="The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of "
                          "the following values: [unsorted, queryname, coordinate, duplicate, unknown]"),
            *Gatk4SortSamBase.additional_args
        ]

    def outputs(self):
        return [
            ToolOutput("output", BamPair(), glob="$(inputs.outputFilename)")
        ]

    additional_args = [
        ToolInput("argumentsFile", Array(File(), optional=True), prefix="--arguments_file", position=10,
                  doc="read one or more arguments files and add them to the command line"),
        ToolInput("compressionLevel", Int(optional=True), prefix="--COMPRESSION_LEVEL", position=11,
                  doc="Compression level for all compressed files created (e.g. BAM and GELI)."),
        ToolInput("createIndex", Boolean(optional=True), prefix="--CREATE_INDEX", position=11,
                  doc="Whether to create a BAM index when writing a coordinate-sorted BAM file."),
        ToolInput("createMd5File", Boolean(optional=True), prefix="--CREATE_MD5_FILE", position=11,
                  doc="Whether to create an MD5 digest for any BAM or FASTQ files created."),
        ToolInput("maxRecordsInRam", Int(optional=True), prefix="--MAX_RECORDS_IN_RAM", position=11,
                  doc="When writing SAM files that need to be sorted, this will specify the number of "
                      "records stored in RAM before spilling to disk. Increasing this number reduces "
                      "the number of file handles needed to sort a SAM file, and increases the amount of RAM needed."),
        ToolInput("quiet", Boolean(optional=True), prefix="--QUIET", position=11,
                  doc="Whether to suppress job-summary info on System.err."),
        ToolInput("reference", FastaWithDict(optional=True), prefix="--reference", position=11,
                  doc="Reference sequence file."),
        ToolInput("tmpDir", Directory(optional=True), prefix="--TMP_DIR", position=11,
                  doc="Undocumented option"),
        ToolInput("useJdkDeflater", Boolean(optional=True), prefix="--use_jdk_deflater", position=11,
                  doc="Whether to use the JdkDeflater (as opposed to IntelDeflater)"),
        ToolInput("useJdkInflater", Boolean(optional=True), prefix="--use_jdk_inflater", position=11,
                  doc="Whether to use the JdkInflater (as opposed to IntelInflater)"),
        ToolInput("validationStringency", String(optional=True), prefix="--VALIDATION_STRINGENCY", position=11,
                  doc="Validation stringency for all SAM files read by this program. Setting stringency to SILENT "
                      "can improve performance when processing a BAM file in which variable-length data "
                      "(read, qualities, tags) do not otherwise need to be decoded."
                      "The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), "
                      "which can have one of the following values: [STRICT, LENIENT, SILENT]"),
        ToolInput("verbosity", String(optional=True), prefix="--verbosity", position=11,
                  doc="The --verbosity argument is an enumerated type (LogLevel), which can have "
                      "one of the following values: [ERROR, WARNING, INFO, DEBUG]")
    ]

