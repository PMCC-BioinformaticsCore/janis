from Pipeline import Array, String, Int, File, CommandTool, ToolOutput, \
    ToolInput, ToolArgument, Double, Boolean, Filename
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.bai import Bai
from Pipeline.bioinformatics.data_types.bam import Bam


class PicardMarkDup(CommandTool):
    inputBam = ToolInput("inputBam", Bam(), position=4,  # type should be Array(Bam(), optional=True)
                         prefix="INPUT=", separate_value_from_prefix=False,
                         doc="One or more input SAM or BAM files to analyze. Must be coordinate sorted. "
                             "Default value null. This option may be specified 0 or more times")

    outputFilename = ToolInput("outputFilename", Filename(extension=".bam"), position=5, prefix="OUTPUT=",
                               separate_value_from_prefix=False,
                               doc="The output file to write marked records to Required")

    output = ToolOutput("output", Bam(), glob='$(inputs.outputFilename)')  # o09 | o12
    metric = ToolOutput("metric", File(), glob='$(inputs.metricsFile)')  # o10 | o13
    index = ToolOutput("index", Bai(), glob='$(inputs.outputFilename.replace(".bam", ".bai"))')

    outputPair = ToolOutput("outputPair", BamPair(), glob='$(inputs.outputFilename)')

    @staticmethod
    def tool():
        return "PicardMarkDups"

    @staticmethod
    def base_command():
        return ['java']

    @staticmethod
    def docker():
        return "biocontainers/picard:v2.3.0_cv3"

    # @staticmethod
    # def environment_variables():
    #     return {
    #         "PATH": "/usr/local/bin/:/usr/bin:/bin:/opt/conda/bin"
    #     }

    @staticmethod
    def doc():
        return "picard-BuildBamIndex.cwl is developed for CWL consortium" \
               "Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules. " \
               "All records are then written to the output file with the duplicate records flagged"

    def arguments(self):
        return [
            ToolArgument("/opt/conda/share/picard-2.3.0-0/picard.jar", position=2, prefix="-jar"),
            ToolArgument("MarkDuplicates", position=3)
        ]

    comment = ToolInput("comment", Array(File(), optional=True), position=17,
                        doc="Comment(s) to include in the output files header. Default value null. "
                            "This option may be specified 0 or more times")
    groupCommandName = ToolInput("groupCommandName", String(optional=True), position=16, prefix="PROGRAM_GROUP_NAME=",
                                 separate_value_from_prefix=False,
                                 doc="Value of PN tag of PG record to be created. Default value MarkDuplicates. "
                                     "This option can be set to 'null' to clear the default value")
    validation_stringency = ToolInput("validation_stringency", String(optional=True), position=3,
                                      prefix="VALIDATION_STRINGENCY=", separate_value_from_prefix=False)
    groupVersion = ToolInput("groupVersion", String(optional=True), position=14, prefix="PROGRAM_GROUP_VERSION=",
                             separate_value_from_prefix=False,
                             doc="Value of VN tag of PG record to be created. If not specified, "
                                 "the version will be detected automatically. Default value null")
    readSorted = ToolInput("readSorted", Boolean(optional=True), position=22, prefix="ASSUME_SORTED=",
                           separate_value_from_prefix=False,
                           doc="If true, assume that the input file is coordinate sorted even if the header says otherwise. "
                               "Default value false. This option can be set to 'null' to clear the default value. "
                               "Possible values {true, false}")
    readOneBarcodeTag = ToolInput("readOneBarcodeTag", String(optional=True), position=11,
                                  prefix="READ_ONE_BARCODE_TAG=", separate_value_from_prefix=False,
                                  doc="Read one barcode SAM tag (ex. BX for 10X Genomics) Default value null")
    metricsFile = ToolInput("metricsFile", Filename(), default="metricsFile-Tumor-markDuplicates",
                            position=6, prefix="METRICS_FILE=", separate_value_from_prefix=False,
                            doc="File to write duplication metrics to Required")
    regularExpression = ToolInput("regularExpression", String(optional=True), position=18, prefix="READ_NAME_REGEX=",
                                  separate_value_from_prefix=False,
                                  doc="Regular expression that can be used to parse read names in the incoming SAM file."
                                      "Read names are parsed to extract three variables tile/region, x coordinate and "
                                      "y coordinate. These values are used to estimate the rate of optical duplication "
                                      "in order to give a more accurate estimated library size. Set this option to null "
                                      "to disable optical duplicate detection. The regular expression should contain "
                                      "three capture groups for the three variables, in order. It must match the entire "
                                      "read name. Note that if the default regex is specified, a regex match is not "
                                      "actually done, but instead the read name is split on colon character. "
                                      "For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, "
                                      "x and y values. For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements "
                                      "are assumed to be tile, x and y values. "
                                      "Default value [a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*. "
                                      "This option can be set to 'null' to clear the default value")
    pixelDistance = ToolInput("pixelDistance", Int(optional=True), position=19,
                              prefix="OPTICAL_DUPLICATE_PIXEL_DISTANCE=", separate_value_from_prefix=False,
                              doc="The maximum offset between two duplicte clusters in order to consider "
                                  "them optical duplicates. This should usually be set to some fairly small number "
                                  "(e.g. 5-10 pixels) unless using later versions of the Illumina pipeline that multiply "
                                  "pixel values by 10, in which case 50-100 is more normal. Default value 100. "
                                  "This option can be set to 'null' to clear the default value")
    picard_markdup_tmpdir = ToolInput("picard_markdup_tmpdir", String(optional=True), position=21, prefix="TMP_DIR=",
                                      separate_value_from_prefix=False,
                                      doc="Default value null. This option may be specified 0 or more times.")
    java_arg = ToolInput("java_arg", String(optional=True), position=1, default="-Xmx4g")
    recordId = ToolInput("recordId", String(optional=True), position=13, prefix="PROGRAM_RECORD_ID=",
                         separate_value_from_prefix=False,
                         doc="The program record ID for the @PG record(s) created by this program. Set to null to "
                             "disable PG record creation. This string may have a suffix appended to avoid collision "
                             "with other program record IDs. Default value MarkDuplicates. This option can be set to "
                             "'null' to clear the default value")
    maxRecordsInRam = ToolInput("maxRecordsInRam", Int(optional=True), position=4, prefix="MAX_RECORDS_IN_RAM=",
                                separate_value_from_prefix=False)
    maxFileHandles = ToolInput("maxFileHandles", Int(optional=True), position=8,
                               prefix="MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=", separate_value_from_prefix=False,
                               doc="Maximum number of file handles to keep open when spilling read ends to disk. "
                                   "Set this number a little lower than the per-process maximum number of file "
                                   "that may be open. This number can be found by executing the 'ulimit -n' command "
                                   "on a Unix system. Default value 8000. This option can be set to 'null' to "
                                   "clear the default value")
    sortRatio = ToolInput("sortRatio", Double(optional=True), position=9, prefix="SORTING_COLLECTION_SIZE_RATIO=",
                          separate_value_from_prefix=False,
                          doc="This number, plus the maximum RAM available to the JVM, determine the memory footprint "
                              "used by some of the sorting collections. If you are running out of memory, try reducing "
                              "this number. Default value 0.25. This option can be set to 'null' to clear the default value")
    groupCommandLine = ToolInput("groupCommandLine", String(optional=True), position=15,
                                 prefix="PROGRAM_GROUP_COMMAND_LINE=", separate_value_from_prefix=False,
                                 doc="Value of CL tag of PG record to be created. If not supplied the command line "
                                     "will be detected automatically. Default value null")
    removeDuplicates = ToolInput("removeDuplicates", Boolean(optional=True), position=7, prefix="REMOVE_DUPLICATES=",
                                 separate_value_from_prefix=False,
                                 doc="If true do not write duplicates to the output file instead of writing them with "
                                     "appropriate flags set. Default value false. This option can be set to 'null' to "
                                     "clear the default value. Possible values {true, false}")
    createIndex = ToolInput("createIndex", String(optional=True), position=20, prefix="CREATE_INDEX=",
                            separate_value_from_prefix=False,
                            doc="Whether to create a BAM index when writing a coordinate-sorted BAM file. "
                                "Default value false. This option can be set to 'null' to clear the default value. "
                                "Possible values {true, false}",
                            default="true")
    readTwoBarcodeTag = ToolInput("readTwoBarcodeTag", String(optional=True), position=12,
                                  prefix="READ_TWO_BARCODE_TAG=", separate_value_from_prefix=False,
                                  doc="Read two barcode SAM tag (ex. BX for 10X Genomics) Default value null")
    barcodeTag = ToolInput("barcodeTag", String(optional=True), position=10, prefix="BARCODE_TAG=",
                           separate_value_from_prefix=False,
                           doc="Barcode SAM tag (ex. BC for 10X Genomics) Default value null")


if __name__ == "__main__":
    print(PicardMarkDup().help())
