from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline import File, String, Array, CommandTool, ToolOutput, ToolInput, ToolArgument, Int, Boolean, Double, Filename
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.vcfidx import VcfIdx


class GatkRecalibrator(CommandTool):
    inputBase = ToolInput("inputBase", BamPair(), position=6, prefix="-I",
                          doc="bam file produced after indelRealigner")

    reference = ToolInput("reference", FastaWithDict(), position=5, prefix="-R")

    known = ToolInput("known", Array(VcfIdx(), optional=True), prefix="--knownSites", position=28,
                      doc="Any number of VCF files representing known SNPs and/or indels. "
                          "Could be e.g. dbSNP and/or official 1000 Genomes indel calls. "
                          "SNPs in these files will be ignored unless the --mismatchFraction argument is used.")
    outputFile = ToolInput("outputFile", Filename(extension=".grp"), position=8, prefix="-o",
                                            doc="name of the output file from baseRecalibrator")

    output = ToolOutput("output", File(), glob='$(inputs.outputFile)')

    @staticmethod
    def tool():
        return "GatkRecalibrator"

    @staticmethod
    def base_command():
        return ['java']

    @staticmethod
    def docker():
        return "broadinstitute/gatk3:3.7-0"

    @staticmethod
    def doc():
        return "GATK-BaseRecalibrator.cwl is developed for CWL consortiumIt generate base recalibration table to " \
               "compensate for systematic errors in basecalling confidences. " \
               "Usage: java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R reference.fasta -I my_reads.bam " \
               "-knownSites latest_dbsnp.vcf -o recal_data.table."

    def arguments(self):
        return [
            ToolArgument("./test/test-files", position=2, prefix="-Djava.io.tmpdir=", separate_value_from_prefix=False),
            ToolArgument("/usr/GenomeAnalysisTK.jar", position=3, prefix="-jar"),
            ToolArgument("BaseRecalibrator", position=4, prefix="-T"),
            ToolArgument("--filter_bases_not_stored", position=30)
        ]

    inputIndex = ToolInput("inputIndex", File(optional=True))
    deletions_default_quality = ToolInput("deletions_default_quality", Int(optional=True), position=24,
                                          prefix="--deletions_default_quality",
                                          doc="default quality for the base deletions covariate")
    binary_tag_name = ToolInput("binary_tag_name", String(optional=True), position=27, prefix="--binary_tag_name",
                                doc="the binary tag covariate name if using it")
    no_standard_covs = ToolInput("no_standard_covs", Boolean(optional=True), position=15, prefix="--no_standard_covs",
                                 doc="Do not use the standard set of covariates, "
                                     "but rather just the ones listed using the -cov argument")
    solid_nocall_strategy = ToolInput("solid_nocall_strategy", String(optional=True), position=11,
                                      prefix="--solid_nocall_strategy",
                                      doc="Defines the behavior of the recalibrator when it encounters "
                                          "no calls in the color space. "
                                          "Options = THROW_EXCEPTION, LEAVE_READ_UNRECALIBRATED, or PURGE_READ")
    low_quality_tail = ToolInput("low_quality_tail", Int(optional=True), position=20, prefix="--low_quality_tail",
                                 doc="minimum quality for the bases in the tail of the reads to be considered")
    java_arg = ToolInput("java_arg", String(optional=True), position=1, default="-Xmx4g")
    out = ToolInput("out", File(optional=True), position=14, prefix="--out",
                    doc="The output recalibration table file to create")
    quantizing_levels = ToolInput("quantizing_levels", Boolean(optional=True), position=13,
                                  prefix="--quantizing_levels",
                                  doc="Sort the rows in the tables of reports. Whether GATK report tables should "
                                      "have rows in sorted order, starting from leftmost column")
    bqsrBAQGapOpenPenalty = ToolInput("bqsrBAQGapOpenPenalty", Double(optional=True), position=26,
                                      prefix="--bqsrBAQGapOpenPenalty",
                                      doc="BQSR BAQ gap open penalty (Phred Scaled). Default value is 40. 30 is "
                                          "perhaps better for whole genome call sets")
    mismatches_context_size = ToolInput("mismatches_context_size", Int(optional=True), position=17,
                                        prefix="--mismatches_context_size",
                                        doc="Size of the k-mer context to be used for base mismatches")
    maximum_cycle_value = ToolInput("maximum_cycle_value", Int(optional=True), position=18,
                                    prefix="--maximum_cycle_value",
                                    doc="The maximum cycle value permitted for the Cycle covariate")
    run_without_dbsnp_potentially_ruining_quality = ToolInput("run_without_dbsnp_potentially_ruining_quality",
                                                              Boolean(optional=True), position=12,
                                                              prefix="--run_without_dbsnp_potentially_ruining_quality",
                                                              doc="If specified, allows the recalibrator to be used "
                                                                  "without a dbsnp rod. Very unsafe and for expert users only.")
    lowMemoryMode = ToolInput("lowMemoryMode", Boolean(optional=True), position=19, prefix="--lowMemoryMode",
                              doc="Reduce memory usage in multi-threaded code at the expense of threading efficiency")
    solid_recal_mode = ToolInput("solid_recal_mode", String(optional=True), position=10, prefix="--solid_recal_mode",
                                 doc="How should we recalibrate solid bases in which the reference was inserted? "
                                     "Options = DO_NOTHING, SET_Q_ZERO, SET_Q_ZERO_BASE_N, or REMOVE_REF_BIAS")
    insertions_default_quality = ToolInput("insertions_default_quality", Int(optional=True), position=22,
                                           prefix="--insertions_default_quality",
                                           doc="default quality for the base insertions covariate")
    sort_by_all_columns = ToolInput("sort_by_all_columns", Boolean(optional=True), position=9,
                                    prefix="--sort_by_all_columns",
                                    doc="Sort the rows in the tables of reports. Whether GATK report tables "
                                        "should have rows in sorted order, starting from leftmost column")
    list = ToolInput("list", Boolean(optional=True), position=21, prefix="--list",
                     doc="List the available covariates and exit")
    indels_context_size = ToolInput("indels_context_size", Int(optional=True), position=23,
                                    prefix="--indels_context_size",
                                    doc="Size of the k-mer context to be used for base insertions and deletions")
    mismatches_default_quality = ToolInput("mismatches_default_quality", Int(optional=True), position=16,
                                           prefix="--mismatches_default_quality",
                                           doc="default quality for the base mismatches covariate")
    downsamplingType = ToolInput("downsamplingType", String(optional=True), position=27, prefix="--downsampling_type",
                                 default="none")
    bedFile = ToolInput("bedFile", Bed(optional=True), position=28, prefix="-L")
    threads = ToolInput("threads", Int(optional=True), position=26, prefix="-nct", default=4)


if __name__ == "__main__":
    print(GatkRecalibrator().help())
