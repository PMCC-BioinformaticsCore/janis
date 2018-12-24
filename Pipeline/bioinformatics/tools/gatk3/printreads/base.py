from abc import ABC

from Pipeline import ToolInput, File, ToolOutput, ToolArgument, Array, String, Boolean, Int
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.tools.gatk3.gatk3base import Gatk3Base


class Gatk3PrintReadsBase(Gatk3Base, ABC):


    @staticmethod
    def tool():
        return "GatkPrintReads"

    def inputs(self):
        return [
            *super(Gatk3PrintReadsBase, self).inputs(),
            ToolInput("inputBam", BamPair(), position=6, prefix="-I",
                      doc="BAM/SAM/CRAM file containing reads"),

            ToolInput("reference", FastaWithDict(), position=5, prefix="-R"),
            ToolInput("input_baseRecalibrator", File(), position=7, prefix="-BQSR",
                      doc="the recalibration table produced by BaseRecalibration"),

            ToolInput("bedFile", Bed(), position=15, prefix="-L"),
            *Gatk3PrintReadsBase.additional_args
        ]

    def outputs(self):
        return [
            # ToolOutput("output_printReads", Bam(), glob='$(inputs.outputfile_printReads)'),
            # ToolOutput("printreads_index_output", Bai(),
            #            glob='$(inputs.outputfile_printReads.replace(".bam", ".bai"))'),
            ToolOutput("pairedOutput", BamPair(), glob='$(inputs.outputFilename)', doc="Write output to this file")

        ]

    @staticmethod
    def doc():
        return """
Write reads from SAM format file (SAM/BAM/CRAM) that pass criteria to a new file.
A common use case is to subset reads by genomic interval using the -L argument. 
Note when applying genomic intervals, the tool is literal and does not retain mates 
of paired-end reads outside of the interval, if any. Data with missing mates will fail 
ValidateSamFile validation with MATE_NOT_FOUND, but certain tools may still analyze the data. 
If needed, to rescue such mates, use either FilterSamReads or ExtractOriginalAlignmentRecordsByNameSpark.

By default, PrintReads applies the WellformedReadFilter at the engine level. What this means is that 
the tool does not print reads that fail the WellformedReadFilter filter. You can similarly apply 
other engine-level filters to remove specific types of reads with the --read-filter argument. 
See documentation category 'Read Filters' for a list of available filters. 
To keep reads that do not pass the WellformedReadFilter, either disable the filter 
with --disable-read-filter or disable all default filters with --disable-tool-default-read-filters.

The reference is strictly required when handling CRAM files.

Documentation: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_PrintReads.php
        """.strip()
        # return "GATK-RealignTargetCreator.cwl is developed for CWL consortiumPrints all reads that have a mapping " \
        #        "quality above zero  Usage: java -Xmx4g -jar GenomeAnalysisTK.jar -T PrintReads -R reference.fasta " \
        #        "-I input1.bam -I input2.bam -o output.bam --read_filter MappingQualityZero"

    def arguments(self):
        return [
            ToolArgument("./test/test-files", position=2, prefix="-Djava.io.tmpdir=", separate_value_from_prefix=False),
            ToolArgument("PrintReads", position=4, prefix="-T"),
            ToolArgument("--filter_bases_not_stored", position=20)
        ]

    additional_args = [
        ToolInput("sample_file", Array(File(), optional=True), position=11),
        ToolInput("platform", String(optional=True), position=13, prefix="--platform",
                             doc="Exclude all reads with this platform from the output"),
        ToolInput("number", String(optional=True), position=13, prefix="--number",
                           doc="Exclude all reads with this platform from the output"),
        ToolInput("simplify", Boolean(optional=True), position=9, prefix="--simplify",
                             doc="Erase all extra attributes in the read but keep the read group information"),
        ToolInput("readGroup", String(optional=True), position=12, prefix="--readGroup",
                              doc="Exclude all reads with this read group from the output"),
        ToolInput("sample_name", Array(String(), optional=True), position=10,
                                doc="Sample name to be included in the analysis. Can be specified multiple times."),
        ToolInput("outputFilename", String(optional=True), position=8, prefix="-o",
                                          doc="name of the output file from indelRealigner"),
        ToolInput("java_arg", String(), position=1, default="-Xmx4g"),
        ToolInput("threads", Int(), position=14, prefix="-nct", default=4),
        ToolInput("downsamplingType", String(optional=True), position=16, prefix="--downsampling_type",
                                     default="none")
    ]


if __name__ == "__main__":
    print(Gatk3PrintReadsBase().help())
