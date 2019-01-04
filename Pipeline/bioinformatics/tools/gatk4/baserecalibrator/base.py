from abc import ABC

from Pipeline import ToolArgument, ToolInput, ToolOutput, Filename, Array, File, String, Int, Boolean, Double, Directory
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.vcf import VcfIdx
from Pipeline.bioinformatics.tools.gatk4.gatk4toolbase import Gatk4ToolBase
from Pipeline.unix.data_types.tsv import Tsv


class Gatk4BaseRecalibratorBase(Gatk4ToolBase, ABC):

    @classmethod
    def gatk_command(cls):
        return "BaseRecalibrator"

    @staticmethod
    def tool():
        return "Gatk4BaseRecalibrator"

    def inputs(self):
        return [
            *super(Gatk4BaseRecalibratorBase, self).inputs(),
            *Gatk4BaseRecalibratorBase.additional_args,

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

            ToolInput("outputFilename", Filename(extension=".table"), position=8, prefix="-O",
                      doc="**The output recalibration table filename to create.** "
                          "After the header, data records occur one per line until the end of the file. The first "
                          "several items on a line are the values of the individual covariates and will change "
                          "depending on which covariates were specified at runtime. The last three items are the "
                          "data- that is, number of observations for this combination of covariates, number of "
                          "reference mismatches, and the raw empirical quality score calculated by phred-scaling "
                          "the mismatch rate. Use '/dev/stdout' to print to standard out.")
        ]

    def outputs(self):
        return [
            ToolOutput("output", Tsv(), glob="$(inputs.outputFilename)")
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
            *super(Gatk4BaseRecalibratorBase, self).arguments(),
            ToolArgument("--filter_bases_not_stored", position=30)
        ]

    additional_args = [
        ToolInput("tmpDir", Directory(optional=True), prefix="--tmp-dir", doc="Temp directory to use.")
    ]


if __name__ == "__main__":
    print(Gatk4BaseRecalibratorBase().help())