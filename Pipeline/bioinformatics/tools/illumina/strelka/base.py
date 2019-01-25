from abc import ABC, abstractmethod
from typing import List

from Pipeline import CommandTool, ToolOutput, ToolInput, ToolArgument, Boolean, String, File
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.fasta import FastaWithDict
from Pipeline.bioinformatics.data_types.vcf import TabixIdx
from Pipeline.unix.data_types.tsv import Tsv
from Pipeline.utils.metadata import ToolMetadata


class StrelkaBase(CommandTool, ABC):
    @staticmethod
    def tool():
        return "strelka-germline"

    @staticmethod
    def base_command():
        return None

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("bam", BamPair(), prefix="--bam", position=1, shell_quote=False,
                      doc="Sample BAM or CRAM file. May be specified more than once, multiple inputs will be treated "
                          "as each BAM file representing a different sample. [required] (no default)"),
            ToolInput("reference", FastaWithDict(), prefix="--referenceFasta", position=1, shell_quote=False,
                      doc="samtools-indexed reference fasta file [required]"),
            ToolInput("relativeStrelkaDirectory", String(optional=True, default="strelka_dir"), prefix="--runDir",
                      position=1, shell_quote=False,
                      doc="Name of directory to be created where all workflow scripts and output will be written. "
                          "Each analysis requires a separate directory."),

            ToolInput("version", Boolean(optional=True), prefix="--version", position=3, shell_quote=False,
                      doc="show program's version number and exit"),
            ToolInput("help", Boolean(optional=True), prefix="--help", position=3, shell_quote=False,
                      doc="(-h) show this help message and exit"),
            ToolInput("mode", String(optional=True, default="local"), prefix="--mode", position=3, shell_quote=False,
                      doc="(-m MODE)  select run mode (local|sge)"),
            ToolInput("queue", String(optional=True), prefix="--queue", position=3, shell_quote=False,
                      doc="(-q QUEUE) specify scheduler queue name"),
            ToolInput("jobs", String(optional=True), prefix="--jobs", position=3, shell_quote=False,
                      doc=" (-j JOBS)  number of jobs, must be an integer or 'unlimited' "
                          "(default: Estimate total cores on this node for local mode, 128 for sge mode)"),
            ToolInput("memGb", String(optional=True), prefix="--memGb", position=3, shell_quote=False,
                      doc=" (-g MEMGB) gigabytes of memory available to run workflow "
                          "-- only meaningful in local mode, must be an integer (default: Estimate the total "
                          "memory for this node for local mode, 'unlimited' for sge mode)"),
            ToolInput("dryRun", Boolean(optional=True), prefix="--dryRun", position=3, shell_quote=False,
                      doc="dryRun (-d,) workflow code without actually running command-tasks"),
            ToolInput("quiet", Boolean(optional=True), prefix="--quiet", position=3, shell_quote=False,
                      doc="Don't write any log output to stderr "
                          "(but still write to workspace/pyflow.data/logs/pyflow_log.txt)"),
            ToolInput("mailTo", String(optional=True), prefix="--mailTo", position=3, shell_quote=False,
                      doc="(-e) send email notification of job completion status to this address "
                          "(may be provided multiple times for more than one email address)"),
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            # ToolOutput("directory", Directory(), glob="$(inputs.relativeStrelkaDirectory)"),
            ToolOutput("configPickle", File(),
                       glob="$(inputs.relativeStrelkaDirectory + '/runWorkflow.py.config.pickle')"),
            ToolOutput("script", File(), glob="$(inputs.relativeStrelkaDirectory + '/runWorkflow.py')"),
            ToolOutput("stats", Tsv(), glob="$(inputs.relativeStrelkaDirectory + '/results/stats/runStats.tsv')",
                       doc="A tab-delimited report of various internal statistics from the variant calling process: "
                           "Runtime information accumulated for each genome segment, excluding auxiliary steps such "
                           "as BAM indexing and vcf merging. Indel candidacy statistics"),
            ToolOutput("variants", TabixIdx(),
                       glob="$(inputs.relativeStrelkaDirectory + '/results/variants/variants.vcf.gz')",
                       doc="Primary variant inferences are provided as a series of VCF 4.1 files"),
            ToolOutput("genome", TabixIdx(),
                       glob="$(inputs.relativeStrelkaDirectory + '/results/variants/genome.vcf.gz')"),
        ]

    def arguments(self) -> List[ToolArgument]:
        return [
            ToolArgument("configureStrelkaGermlineWorkflow.py", position=0, shell_quote=False),
            ToolArgument(";$(inputs.relativeStrelkaDirectory + '/runWorkflow.py')", position=2, shell_quote=False)

        ]

    @staticmethod
    def requirements():
        from cwlgen.cwlgen import ShellCommandRequirement
        return [ShellCommandRequirement()]

    @staticmethod
    @abstractmethod
    def docker():
        raise Exception("Strelka version must override docker command")

    def friendly_name(self):
        return "Strelka (Germline)"

    def metadata(self):
        from datetime import date
        return ToolMetadata(
            creator="Michael Franklin",
            maintainer="Michael Franklin",
            maintainer_email="michael.franklin@petermac.org",
            date_created=date(2018, 12, 24),
            date_updated=date(2019, 1, 24),
            institution="Illumina",
            doi=None,
            citation=None, # find citation
            keywords=["broad", "igvtools", "index"],
            documentation_url="https://github.com/Illumina/strelka",
            documentation="""
Strelka2 is a fast and accurate small variant caller optimized for analysis of germline variation 
in small cohorts and somatic variation in tumor/normal sample pairs. The germline caller employs 
an efficient tiered haplotype model to improve accuracy and provide read-backed phasing, adaptively 
selecting between assembly and a faster alignment-based haplotyping approach at each variant locus. 
The germline caller also analyzes input sequencing data using a mixture-model indel error estimation 
method to improve robustness to indel noise. The somatic calling model improves on the original 
Strelka method for liquid and late-stage tumor analysis by accounting for possible tumor cell 
contamination in the normal sample. A final empirical variant re-scoring step using random forest 
models trained on various call quality features has been added to both callers to further improve precision.

Compared with submissions to the recent PrecisonFDA Consistency and Truth challenges, the average 
indel F-score for Strelka2 running in its default configuration is 3.1% and 0.08% higher, respectively, 
than the best challenge submissions. Runtime on a 28-core server is ~40 minutes for 40x WGS germline 
analysis and ~3 hours for a 110x/40x WGS tumor-normal somatic analysis

Strelka accepts input read mappings from BAM or CRAM files, and optionally candidate and/or forced-call 
alleles from VCF. It reports all small variant predictions in VCF 4.1 format. Germline variant 
reporting uses the gVCF conventions to represent both variant and reference call confidence. 
For best somatic indel performance, Strelka is designed to be run with the Manta structural variant 
and indel caller, which provides additional indel candidates up to a given maxiumum indel size 
(49 by default). By design, Manta and Strelka run together with default settings provide complete 
coverage over all indel sizes (in additional to SVs and SNVs). 

See the user guide for a full description of capabilities and limitations""".strip()
        )
