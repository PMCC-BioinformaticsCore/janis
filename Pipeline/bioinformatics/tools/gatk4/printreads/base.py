from abc import ABC
from typing import List

from Pipeline import ToolOutput, ToolInput, Filename
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.tools.gatk4.gatk4toolbase import Gatk4ToolBase


class Gatk4PrintReadsBase(Gatk4ToolBase, ABC):

    @classmethod
    def gatk_command(cls):
        return "PrintReads"

    @staticmethod
    def tool():
        return "Gatk4PrintReads"

    def inputs(self):
        return [
            ToolInput("input", Bam()),
            ToolInput("outputFilename", Filename())
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("output", BamPair(), glob="$(inputs.outputFilename)")
        ]


    @staticmethod
    def doc():
        return """
        Overview
Write reads from SAM format file (SAM/BAM/CRAM) that pass criteria to a new file.
A common use case is to subset reads by genomic interval using the -L argument. Note when applying genomic intervals, the tool is literal and does not retain mates of paired-end reads outside of the interval, if any. Data with missing mates will fail ValidateSamFile validation with MATE_NOT_FOUND, but certain tools may still analyze the data. If needed, to rescue such mates, use either FilterSamReads or ExtractOriginalAlignmentRecordsByNameSpark.

By default, PrintReads applies the WellformedReadFilter at the engine level. What this means is that the tool does not print reads that fail the WellformedReadFilter filter. You can similarly apply other engine-level filters to remove specific types of reads with the --read-filter argument. See documentation category 'Read Filters' for a list of available filters. To keep reads that do not pass the WellformedReadFilter, either disable the filter with --disable-read-filter or disable all default filters with --disable-tool-default-read-filters.

The reference is strictly required when handling CRAM files.

Documentation: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_PrintReads.php
        
        """