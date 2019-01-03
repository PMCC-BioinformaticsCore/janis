from abc import ABC

from Pipeline import ToolInput, Filename, ToolOutput, File, String, Directory
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.tools.gatk4.gatk4toolbase import Gatk4ToolBase
from Pipeline.unix.data_types.tsv import Tsv


class Gatk4ApplyBqsrBase(Gatk4ToolBase, ABC):
    @staticmethod
    def tool():
        return "GATK4ApplyBQSR"

    def inputs(self):
        return [
            ToolInput("input", BamPair(), prefix="-I", doc="The SAM/BAM/CRAM file containing reads.", position=10),
            ToolInput("reference", FastaWithDict(), prefix="-R", doc="Reference sequence"),
            ToolInput("outputFilename", Filename(extension="recal.bam"), prefix="-O", doc="Write output to this file"),
            ToolInput("recalFile", Tsv(optional=True), prefix="--bqsr-recal-file",
                      doc="Input recalibration table for BQSR"),

            *super(Gatk4ApplyBqsrBase, self).inputs(),
            *self.additional_args
        ]

    def outputs(self):
        return [
            ToolOutput("output", BamPair(), glob="$(inputs.outputFilename)")
        ]

    @staticmethod
    def doc():
        return """
    Apply base quality score recalibration: This tool performs the second pass in a two-stage 
    process called Base Quality Score Recalibration (BQSR). Specifically, it recalibrates the 
    base qualities of the input reads based on the recalibration table produced by the 
    BaseRecalibrator tool, and outputs a recalibrated BAM or CRAM file.

    Summary of the BQSR procedure: The goal of this procedure is to correct for systematic bias 
    that affect the assignment of base quality scores by the sequencer. The first pass consists 
    of calculating error empirically and finding patterns in how error varies with basecall 
    features over all bases. The relevant observations are written to a recalibration table. 
    The second pass consists of applying numerical corrections to each individual basecall 
    based on the patterns identified in the first step (recorded in the recalibration table) 
    and write out the recalibrated data to a new BAM or CRAM file.
    
    - This tool replaces the use of PrintReads for the application of base quality score 
        recalibration as practiced in earlier versions of GATK (2.x and 3.x).
    - You should only run ApplyBQSR with the covariates table created from the input BAM or CRAM file(s).
    - Original qualities can be retained in the output file under the "OQ" tag if desired. 
        See the `--emit-original-quals` argument for details.
    
    Documentation: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php
    """.strip()

    additional_args = [
        # Put more detail in here from documentation
        ToolInput("tmpDir", Directory(optional=True), prefix="--TMP_DIR", position=11,
                  doc="Undocumented option"),
    ]
