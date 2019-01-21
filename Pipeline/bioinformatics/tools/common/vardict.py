from typing import List

from Pipeline import CommandTool, ToolOutput, ToolInput, Filename, File, ToolArgument, Boolean, Float, Int, String
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.fasta import FastaFai
from Pipeline.bioinformatics.data_types.vcf import Vcf
from Pipeline.types.common_data_types import Stdout


class VarDict(CommandTool):
    @staticmethod
    def tool():
        return "SplitMultiAllele"

    @staticmethod
    def base_command():
        return [
            "export VarDict=\"/config/binaries/vardict/1.5.1/bin/VarDict\";",
            "$VarDict"
        ]

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("input", Bed(), position=2),
            ToolInput("outputFilename", Filename(extension=".vcf")),
            ToolInput("indexedBam", String(), prefix="-b", position=1, doc="The indexed BAM file"),

            ToolInput("referenceFasta", FastaFai(), prefix="-G", position=1,
                      doc="The reference fasta. Should be indexed (.fai). "
                          "Defaults to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa"),
            # ToolInput("reference", FastaWithDict()) # eventually
            *VarDict.vardict_inputs,
            *VarDict.var2vcf_inputs
        ]

    def outputs(self):
        return [
            ToolOutput("output", Stdout(Vcf()))
        ]

    def arguments(self):
        return [
            ToolArgument("teststrandbias.R |", position=3, shell_quote=False),
            ToolArgument("var2vcf_valid.pl", position=4, shell_quote=False)
        ]

    @staticmethod
    def requirements():
        from cwlgen.cwlgen import ShellCommandRequirement
        return [ShellCommandRequirement()]

    vardict_inputs = [
        ToolInput("indels3prime", Boolean(optional=True), prefix="-3", position=1,
                  doc="Indicate to move indels to 3-prime if alternative alignment can be achieved."),
        ToolInput("amplicon", Float(optional=True), prefix="-a", position=1,
                  doc="Indicate it's amplicon based calling.  Reads that don't map to the amplicon will be skipped.  "
                      "A read pair is considered belonging  to the amplicon if the edges are less than int bp to "
                      "the amplicon, and overlap fraction is at least float.  Default: 10:0.95"),
        ToolInput("minReads", Int(optional=True), prefix="-B", position=1,
                  doc="The minimum # of reads to determine strand bias, default 2"),
        ToolInput("chromNamesAreNumbers", Boolean(optional=True), prefix="-C", position=1,
                  doc="Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2"),
        ToolInput("chromColumn", Int(optional=True), prefix="-c", position=1, doc="The column for chromosome"),
        ToolInput("debug", Boolean(optional=True), prefix="-D", position=1,
                  doc="Debug mode.  Will print some error messages and append full genotype at the end."),
        ToolInput("splitDelimeter", String(optional=True), prefix="-d", position=1,
                  doc="The delimiter for split region_info, default to tab \"\t\""),
        ToolInput("geneEndCol", Int(optional=True), prefix="-E", position=1,
                  doc="The column for region end, e.g. gene end"),
        ToolInput("segEndCol", Int(optional=True), prefix="-e", position=1,
                  doc="The column for segment ends in the region, e.g. exon ends"),
        ToolInput("filter", String(optional=True), prefix="-F", position=1,
                  doc="The hexical to filter reads using samtools. Default: 0x500 (filter 2nd alignments and "
                      "duplicates). Use -F 0 to turn it off."),
        ToolInput("alleleFreqThreshold", Int(optional=True), prefix="-f", position=1,
                  doc="The threshold for allele frequency, default: 0.05 or 5%"),
        ToolInput("geneNameCol", Int(optional=True), prefix="-g", position=1,
                  doc="The column for gene name, or segment annotation"),
        # ToolInput("help", Boolean(optional=True), prefix="-H", position=1, doc="Print this help page"),
        ToolInput("printHeaderRow", Boolean(optional=True), prefix="-h", position=1,
                  doc="Print a header row describing columns"),
        ToolInput("indelSize", Int(optional=True), prefix="-I", position=1, doc="The indel size.  Default: 120bp"),
        ToolInput("outputSplice", Boolean(optional=True), prefix="-i", position=1, doc="Output splicing read counts"),
        ToolInput("performLocalRealignment", Int(optional=True), prefix="-k", position=1,
                  doc="Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it. "
                      "For Ion or PacBio, 0 is recommended."),
        ToolInput("minMatches", Int(optional=True), prefix="-M", position=1,
                  doc="The minimum matches for a read to be considered. If, after soft-clipping, the matched "
                      "bp is less than INT, then the read is discarded. It's meant for PCR based targeted sequencing "
                      "where there's no insert and the matching is only the primers. Default: 0, or no filtering"),
        ToolInput("maxMismatches", Int(optional=True), prefix="-m", position=1,
                  doc="If set, reads with mismatches more than INT will be filtered and ignored. "
                      "Gaps are not counted as mismatches. Valid only for bowtie2/TopHat or BWA aln "
                      "followed by sampe. BWA mem is calculated as NM - Indels. "
                      "Default: 8, or reads with more than 8 mismatches will not be used."),
        ToolInput("sampleName", String(optional=True), prefix="-N", position=1,
                  doc="The sample name to be used directly.  Will overwrite -n option"),
        ToolInput("regexSampleName", String(optional=True), prefix="-n", position=1,
                  doc="The regular expression to extract sample name from BAM filenames. "
                      "Default to: /([^\/\._]+?)_[^\/]*.bam/"),
        ToolInput("mapq", String(optional=True), prefix="-O", position=1,
                  doc="The reads should have at least mean MapQ to be considered a valid variant. "
                      "Default: no filtering"),
        ToolInput("qratio", Float(optional=True), prefix="-o", position=1,
                  doc="The Qratio of (good_quality_reads)/(bad_quality_reads+0.5). "
                      "The quality is defined by -q option.  Default: 1.5"),
        ToolInput("readPosition", Float(optional=True), prefix="-P", position=1,
                  doc="The read position filter. If the mean variants position is less that specified, "
                      "it's considered false positive.  Default: 5"),
        ToolInput("pileup", Boolean(optional=True), prefix="-p", position=1,
                  doc="Do pileup regardless of the frequency"),
        ToolInput("minMappingQual", Int(optional=True), prefix="-Q", position=1,
                  doc="If set, reads with mapping quality less than INT will be filtered and ignored"),
        ToolInput("phredScore", Int(optional=True), prefix="-q", position=1,
                  doc="The phred score for a base to be considered a good call.  "
                      "Default: 25 (for Illumina) For PGM, set it to ~15, as PGM tends to under estimate base quality."),
        ToolInput("region", String(optional=True), prefix="-R", position=1,
                  doc="The region of interest.  In the format of chr:start-end.  If end is omitted, "
                      "then a single position.  No BED is needed."),
        ToolInput("minVariantReads", Int(optional=True), prefix="-r", position=1,
                  doc="The minimum # of variant reads, default 2"),
        ToolInput("regStartCol", Int(optional=True), prefix="-S", position=1,
                  doc="The column for region start, e.g. gene start"),
        ToolInput("segStartCol", Int(optional=True), prefix="-s", position=1,
                  doc="The column for segment starts in the region, e.g. exon starts"),
        ToolInput("minReadsBeforeTrim", Int(optional=True), prefix="-T", position=1,
                  doc="Trim bases after [INT] bases in the reads"),
        ToolInput("removeDuplicateReads", Boolean(optional=True), prefix="-t", position=1,
                  doc="Indicate to remove duplicated reads.  Only one pair with same start positions will be kept"),
        ToolInput("threads", Int(optional=True), prefix="-th", position=1, doc="Threads count."),
        ToolInput("freq", Int(optional=True), prefix="-V", position=1,
                  doc="The lowest frequency in the normal sample allowed for a putative somatic mutation. "
                      "Defaults to 0.05"),
        ToolInput("vcfFormat", Boolean(optional=True), prefix="-v", position=1, doc="VCF format output"),
        ToolInput("vs", String(optional=True), prefix="-VS", position=1,
                  doc="[STRICT | LENIENT | SILENT] How strict to be when reading a SAM or BAM: "
                      "STRICT   - throw an exception if something looks wrong. "
                      "LENIENT	- Emit warnings but keep going if possible. "
                      "SILENT	- Like LENIENT, only don't emit warning messages. "
                      "Default: LENIENT"),
        ToolInput("bp", Int(optional=True), prefix="-X", position=1,
                  doc="Extension of bp to look for mismatches after insersion or deletion.  "
                      "Default to 3 bp, or only calls when they're within 3 bp."),
        ToolInput("extensionNucleotide", Int(optional=True), prefix="-x", position=1,
                  doc="The number of nucleotide to extend for each segment, default: 0"),
        ToolInput("y", Boolean(optional=True), prefix="-y", position=1, doc="<No content>"),
        ToolInput("downsamplingFraction", Int(optional=True), prefix="-Z", position=1,
                  doc="For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  "
                      "Default: No downsampling.  Use with caution.  "
                      "The downsampling will be random and non-reproducible."),
        ToolInput("zeroBasedCoords", Int(optional=True), prefix="-z", position=1,
                  doc="0/1  Indicate whether coordinates are zero-based, as IGV uses.  "
                      "Default: 1 for BED file or amplicon BED file. Use 0 to turn it off. "
                      "When using the -R option, it's set to 0"),
    ]

    var2vcf_inputs = [
        ToolInput("var2vcfSampleName", String()),
        ToolInput("var2vcfAlleleFreqThreshold", Float()),
    ]


if __name__ == "__main__":
    print(VarDict().help())
