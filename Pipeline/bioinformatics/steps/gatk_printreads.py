from Pipeline.bioinformatics.data_types.bai import Bai
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline import File, String, CommandTool, ToolOutput, ToolInput, ToolArgument, Array, Int, Boolean


# class GatkPrintReads(CommandTool):
#     inputBam_printReads = ToolInput("inputBam_printReads", Bam())
#     input_baseRecalibrator = ToolInput("input_baseRecalibrator", File())
#     reference = ToolInput("reference", Fasta())
#     outputfile_printReads = ToolInput("outputfile_printReads", String())
#     bedFile = ToolInput("bedFile", Bed())
#
#     out = ToolOutput("out", File())
#     out_idx = ToolOutput("out_idx", File())
#
#     @staticmethod
#     def tool():
#         return "gatk-printreads"
#
#     @staticmethod
#     def base_command():
#         return "javac"


class GatkPrintReads(CommandTool):
    reference = ToolInput("reference", FastaWithDict(), position=5, prefix="-R")
    input_baseRecalibrator = ToolInput("input_baseRecalibrator", File(), position=7, prefix="-BQSR",
                                       doc="the recalibration table produced by BaseRecalibration")
    inputBam_printReads = ToolInput("inputBam_printReads", BamPair(), position=6, prefix="-I",
                                    doc="bam file produced after indelRealigner")
    bedFile = ToolInput("bedFile", Bed(), position=15, prefix="-L")

    output_printReads = ToolOutput("output_printReads", Bam(), glob='$(inputs.outputfile_printReads)')
    printreads_index_output = ToolOutput("printreads_index_output", Bai(),
                                         glob='$(inputs.outputfile_printReads.replace(".bam", ".bai"))')
    pairedOutput = ToolOutput("pairedOutput", BamPair(), glob='$(inputs.outputfile_printReads)')

    @staticmethod
    def tool():
        return "GatkPrintReads"

    @staticmethod
    def base_command():
        return ['java']

    @staticmethod
    def docker():
        return "broadinstitute/gatk3:3.7-0"

    @staticmethod
    def doc():
        return "GATK-RealignTargetCreator.cwl is developed for CWL consortiumPrints all reads that have a mapping " \
               "quality above zero  Usage: java -Xmx4g -jar GenomeAnalysisTK.jar -T PrintReads -R reference.fasta " \
               "-I input1.bam -I input2.bam -o output.bam --read_filter MappingQualityZero"

    def arguments(self):
        return [
            ToolArgument("./test/test-files", position=2, prefix="-Djava.io.tmpdir=", separate_value_from_prefix=False),
            ToolArgument("/usr/GenomeAnalysisTK.jar", position=3, prefix="-jar"),
            ToolArgument("PrintReads", position=4, prefix="-T"),
            ToolArgument("--filter_bases_not_stored", position=20)
        ]

    sample_file = ToolInput("sample_file", Array(File(), optional=True), position=11)
    platform = ToolInput("platform", String(optional=True), position=13, prefix="--platform",
                         doc="Exclude all reads with this platform from the output")
    number = ToolInput("number", String(optional=True), position=13, prefix="--number",
                       doc="Exclude all reads with this platform from the output")
    simplify = ToolInput("simplify", Boolean(optional=True), position=9, prefix="--simplify",
                         doc="Erase all extra attributes in the read but keep the read group information")
    readGroup = ToolInput("readGroup", String(optional=True), position=12, prefix="--readGroup",
                          doc="Exclude all reads with this read group from the output")
    sample_name = ToolInput("sample_name", Array(String(), optional=True), position=10,
                            doc="Sample name to be included in the analysis. Can be specified multiple times.")
    outputfile_printReads = ToolInput("outputfile_printReads", String(optional=True), position=8, prefix="-o",
                                      doc="name of the output file from indelRealigner")
    java_arg = ToolInput("java_arg", String(), position=1, default="-Xmx4g")
    threads = ToolInput("threads", Int(), position=14, prefix="-nct", default=4)
    downsamplingType = ToolInput("downsamplingType", String(optional=True), position=16, prefix="--downsampling_type",
                                 default="none")


if __name__ == "__main__":
    print(GatkPrintReads().help())
