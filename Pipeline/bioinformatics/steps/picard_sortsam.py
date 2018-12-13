from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline import File, String, Int, CommandTool, ToolOutput, ToolInput, Boolean, ToolArgument
from Pipeline.types.filename import Filename


class PicardSortSam(CommandTool):
    inputFile_sortSam = ToolInput("inputFile_sortSam", File(), position=4, prefix="INPUT=",
                                      separate_value_from_prefix=False, doc="The BAM or SAM file to sort.")
    outputFileName_sortSam = ToolInput("outputFileName_sortSam", Filename(), position=5, prefix="OUTPUT=",
                                       separate_value_from_prefix=False, doc="The sorted BAM or SAM output file.")

    createIndex = ToolInput("createIndex", Boolean(), default=True, position=8, prefix="CREATE_INDEX=",
                            separate_value_from_prefix=False,
                            doc="Whether to create a BAM index when writing a coordinate-sorted BAM file. "
                                "Default value True.  This option can be set to 'null' to clear the default value. "
                                "Possible values {true, false}")

    soCoordinate = ToolInput("soInput", String(), default="coordinate", position=6, prefix="SORT_ORDER=",
                             separate_value_from_prefix=False,
                             doc="Sort order of output file Required. "
                                 "Possible values {unsorted, queryname, coordinate, duplicate}")

    javaArg = ToolInput("javaArg", String(), default="-Xmx2g", position=1)
    maxRecordsInRam = ToolInput("maxRecordsInRam", Int(optional=True), prefix="MAX_RECORDS_IN_RAM=", position=13,
                                separate_value_from_prefix=False)
    validation_stringency = ToolInput("validation_stringency", String(), default=, prefix="VALIDATION_STRINGENCY", position=10,
                                      separate_value_from_prefix=False)

    tmpdir = ToolInput("tmpdir", Filename(), prefix="TMP_DIR", position=7,
                       doc="This option may be specified 0 or more times.")

    out = ToolOutput("out", Bam(), glob="$(inputs.outputFileName_sortSam)")                                  # Bam file
    indexes = ToolOutput("indexes", File(), glob='$(inputs.outputFileName_sortSam.replace(".bam", ".bai"))') # Bai Index

    @staticmethod
    def tool():
        return "picard-sortsam"

    @staticmethod
    def base_command():
        return "java"

    @staticmethod
    def docker():
        return "biocontainers/picard"

    @staticmethod
    def doc():
        return "picard-SortSam.cwl is developed for CWL consortium. Generates a sorted file."

    def arguments(self):
        return [
            ToolArgument("/opt/conda/share/picard-2.3.0-0/picard.jar", position=2, prefix="-jar"),
            ToolArgument("SortSam", position=3)
        ]


if __name__ == "__main__":
    print(PicardSortSam().help())