from abc import ABC

from Pipeline import ToolArgument, ToolInput, ToolOutput, Filename, Array, File, String, Int, Boolean, Double
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.fasta import FastaWithDict
from Pipeline.bioinformatics.data_types.vcf import VcfIdx
from Pipeline.bioinformatics.tools.gatk3.gatk3toolbase import Gatk3ToolBase


class Gatk3RecalibratorBase(Gatk3ToolBase, ABC):
    output = ToolOutput("output", File(), glob='$(inputs.outputFile)')

    @staticmethod
    def tool():
        return "Gatk3BaseRecalibrator"

    @staticmethod
    def analysis_type():
        return "BaseRecalibrator"

    def inputs(self):
        return [
            *super(Gatk3RecalibratorBase, self).inputs(),
            *Gatk3RecalibratorBase.additional_args,

            ToolInput("input", BamPair(), position=6, prefix="-I", doc="BAM/SAM/CRAM file containing reads"),
            ToolInput("knownSites", Array(VcfIdx()), prefix="--knownSites", position=28,
                      doc="**One or more databases of known polymorphic sites used to exclude "
                          "regions around known polymorphisms from analysis.** "
                          "This algorithm treats every reference mismatch as an indication of error. However, real "
                          "genetic variation is expected to mismatch the reference, so it is critical that a "
                          "database of known polymorphic sites is given to the tool in order to skip over those sites. "
                          "This tool accepts any number of Feature-containing files (VCF, BCF, BED, etc.) for use as "
                          "this database. For users wishing to exclude an interval list of known variation simply "
                          "use -XL my.interval.list to skip over processing those sites. Please note however "
                          "that the statistics reported by the tool will not accurately reflected those sites "
                          "skipped by the -XL argument."),
            ToolInput("reference", FastaWithDict(), position=5, prefix="-R", doc="Reference sequence file"),

            ToolInput("outputFile", Filename(extension=".grp"), position=8, prefix="-O",
                      doc="**The output recalibration table filename to create.** "
                          "After the header, data records occur one per line until the end of the file. The first "
                          "several items on a line are the values of the individual covariates and will change "
                          "depending on which covariates were specified at runtime. The last three items are the "
                          "data- that is, number of observations for this combination of covariates, number of "
                          "reference mismatches, and the raw empirical quality score calculated by phred-scaling "
                          "the mismatch rate. Use '/dev/stdout' to print to standard out.")
        ]

    @staticmethod
    def doc():
        return """
    First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. 
    The default covariates are read group, reported quality score, machine cycle, and nucleotide context.
    
    This walker generates tables based on specified covariates. It does a by-locus traversal operating only at sites 
    that are in the known sites VCF. ExAc, gnomAD, or dbSNP resources can be used as known sites of variation. 
    We assume that all reference mismatches we see are therefore errors and indicative of poor base quality. 
    Since there is a large amount of data one can then calculate an empirical probability of error given the 
    particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a 
    table (of the several covariate values, num observations, num mismatches, empirical quality score).  
    
    Documentation: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php
""".strip()

    def arguments(self):
        return [
            *super(Gatk3RecalibratorBase, self).arguments(),
            ToolArgument("./test/test-files", position=2, prefix="-Djava.io.tmpdir=", separate_value_from_prefix=False),
            ToolArgument("BaseRecalibrator", position=4, prefix="-T"),
            ToolArgument("--filter_bases_not_stored", position=30)
        ]

    additional_args = [
        ToolInput("inputIndex", File(optional=True)),
        ToolInput("deletions_default_quality", Int(optional=True), position=24,
                  prefix="--deletions_default_quality",
                  doc="default quality for the base deletions covariate"),
        ToolInput("binary_tag_name", String(optional=True), position=27, prefix="--binary_tag_name",
                  doc="the binary tag covariate name if using it"),
        ToolInput("no_standard_covs", Boolean(optional=True), position=15, prefix="--no_standard_covs",
                  doc="Do not use the standard set of covariates, "
                      "but rather just the ones listed using the -cov argument"),
        ToolInput("solid_nocall_strategy", String(optional=True), position=11,
                  prefix="--solid_nocall_strategy",
                  doc="Defines the behavior of the recalibrator when it encounters "
                      "no calls in the color space. Options = THROW_EXCEPTION, LEAVE_READ_UNRECALIBRATED, or PURGE_READ"),
        ToolInput("low_quality_tail", Int(optional=True), position=20, prefix="--low_quality_tail",
                  doc="minimum quality for the bases in the tail of the reads to be considered"),
        ToolInput("java_arg", String(optional=True), position=1, default="-Xmx4g"),
        ToolInput("out", File(optional=True), position=14, prefix="--out",
                  doc="The output recalibration table file to create"),
        ToolInput("quantizing_levels", Boolean(optional=True), position=13,
                  prefix="--quantizing_levels",
                  doc="Sort the rows in the tables of reports. Whether GATK report tables should "
                      "have rows in sorted order, starting from leftmost column"),
        ToolInput("bqsrBAQGapOpenPenalty", Double(optional=True), position=26,
                  prefix="--bqsrBAQGapOpenPenalty",
                  doc="BQSR BAQ gap open penalty (Phred Scaled). Default value is 40. 30 is "
                      "perhaps better for whole genome call sets"),
        ToolInput("mismatches_context_size", Int(optional=True), position=17,
                  prefix="--mismatches_context_size",
                  doc="Size of the k-mer context to be used for base mismatches"),
        ToolInput("maximum_cycle_value", Int(optional=True), position=18,
                  prefix="--maximum_cycle_value",
                  doc="The maximum cycle value permitted for the Cycle covariate"),
        ToolInput("run_without_dbsnp_potentially_ruining_quality",
                  Boolean(optional=True), position=12,
                  prefix="--run_without_dbsnp_potentially_ruining_quality",
                  doc="If specified, allows the recalibrator to be used "
                      "without a dbsnp rod. Very unsafe and for expert users only."),
        ToolInput("lowMemoryMode", Boolean(optional=True), position=19, prefix="--lowMemoryMode",
                  doc="Reduce memory usage in multi-threaded code at the expense of threading efficiency"),
        ToolInput("solid_recal_mode", String(optional=True), position=10, prefix="--solid_recal_mode",
                  doc="How should we recalibrate solid bases in which the reference was inserted? "
                      "OPTIONS=DO_NOTHING, SET_Q_ZERO, SET_Q_ZERO_BASE_N, or REMOVE_REF_BIAS"),
        ToolInput("insertions_default_quality", Int(optional=True), position=22,
                  prefix="--insertions_default_quality",
                  doc="default quality for the base insertions covariate"),
        ToolInput("sort_by_all_columns", Boolean(optional=True), position=9,
                  prefix="--sort_by_all_columns",
                  doc="Sort the rows in the tables of reports. Whether GATK report tables "
                      "should have rows in sorted order, starting from leftmost column"),
        ToolInput("list", Boolean(optional=True), position=21, prefix="--list",
                  doc="List the available covariates and exit"),
        ToolInput("indels_context_size", Int(optional=True), position=23,
                  prefix="--indels_context_size",
                  doc="Size of the k-mer context to be used for base insertions and deletions"),
        ToolInput("mismatches_default_quality", Int(optional=True), position=16,
                  prefix="--mismatches_default_quality",
                  doc="default quality for the base mismatches covariate"),
        ToolInput("downsamplingType", String(optional=True), position=27, prefix="--downsampling_type",
                  default="none"),
        ToolInput("bedFile", Bed(optional=True), position=28, prefix="-L"),
        ToolInput("threads", Int(optional=True), position=26, prefix="-nct", default=4)
    ]

if __name__ == "__main__":
    print(Gatk3RecalibratorBase().help())