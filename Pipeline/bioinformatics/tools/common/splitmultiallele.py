from typing import List

from Pipeline import CommandTool, ToolOutput, ToolInput, Filename, File, ToolArgument, Boolean, Float, Int, String
from Pipeline.bioinformatics.data_types.fasta import FastaWithDict
from Pipeline.bioinformatics.data_types.vcf import Vcf
from Pipeline.types.common_data_types import Stdout


class SplitMultiAllele(CommandTool):
    @staticmethod
    def tool():
        return "SplitMultiAllele"

    @staticmethod
    def base_command():
        return "vcfsplitmultiallele.sh"

    @staticmethod
    def docker():
        return "SEE (mfranklin's notes in) DOCUMENTATION"

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("input", Vcf()),
            ToolInput("outputFilename", Filename(extension=".vcf")),
            # ToolInput("reference", FastaWithDict()) # eventually
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("output", Vcf())
        ]

    @staticmethod
    def doc():
        return """
    VcfSplitMultiAllele.sh
    
    Currently stored at: '/researchers/jiaan.yu/WGS_pipeline/VcfSplitMultiAllele.sh’
    
    MFRANKLIN's notes:
        - At the moment, a random shell script is not portable as it exists in a location. I'd be happy
            with putting it in a docker container for when we work out how to run them on the cluster.
            
        - Some of the references are hard-coded, ie: the human reference genome. This won't correctly connect
            when using the execution engine, as I can't reference secondary (index) files and hence the execution
            engine won't correctly localize these.    
        """.strip()



class SplitMultiAllele2(CommandTool):
    @staticmethod
    def tool():
        return "SplitMultiAllele"

    @staticmethod
    def base_command():
        return [
            "export VarDict=\"/config/binaries/vardict/1.5.1/bin/VarDict\";",
            "VarDict -C -v -f 0.05 -N NA12878 -c 1 -S 2 -E 3 -b"
        ]

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("input", Vcf(), position=2),
            ToolInput("outputFilename", Filename(extension=".vcf")),
            # ToolInput("reference", FastaWithDict()) # eventually
        ]

    additional_args = [
        ToolInput("")
    ]

    def outputs(self):
        return [
            ToolOutput("output", Stdout(File()))
        ]

    def arguments(self):
        return [
            ToolArgument("teststrandbias.R |", position=3),
            ToolArgument("var2vcf_Valid.pl")
            ]

    @staticmethod
    def requirements():
        from cwlgen.cwlgen import ShellCommandRequirement
        return [ShellCommandRequirement()]


    vardict_inputs = [
        ToolInput("indels3prime", Boolean(optional=True), prefix="-3", position=1,
                  doc="Indicate to move indels to 3-prime if alternative alignment can be achieved."),
        ToolInput("a", Float(optional=True), prefix="-a", position=1,
                  doc="Indicate it's amplicon based calling.  Reads that don't map to the amplicon will be skipped.  "
                      "A read pair is considered belonging  to the amplicon if the edges are less than int bp to "
                      "the amplicon, and overlap fraction is at least float.  Default: 10:0.95"),
        ToolInput("B", Int(optional=True), prefix="-B", position=1,
                  doc="The minimum # of reads to determine strand bias, default 2"),
        ToolInput("b", String(optional=True), prefix="-b", position=1, doc="The indexed BAM file"),
        ToolInput("C", Boolean(optional=True), prefix="-C", position=1,
                  doc="Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2"),
        ToolInput("c", Int(optional=True), prefix="-c", position=1, doc="The column for chromosome"),
        ToolInput("D", Boolean(optional=True), prefix="-D", position=1,
                  doc="Debug mode.  Will print some error messages and append full genotype at the end."),
        ToolInput("d", String(optional=True), prefix="-d", position=1,
                  doc="The delimiter for split region_info, default to tab \"\t\""),
        ToolInput("E", Int(optional=True), prefix="-E", position=1, doc="The column for region end, e.g. gene end"),
        ToolInput("e", Int(optional=True), prefix="-e", position=1,
                  doc="The column for segment ends in the region, e.g. exon ends"),
        ToolInput("F", String(optional=True), prefix="-F", position=1,
                  doc="The hexical to filter reads using samtools. Default: 0x500 (filter 2nd alignments and "
                      "duplicates). Use -F 0 to turn it off."),
        ToolInput("f", Int(optional=True), prefix="-f", position=1,
                  doc="The threshold for allele frequency, default: 0.05 or 5%"),
        ToolInput("G", String(optional=True), prefix="-G", position=1,
                  doc="The reference fasta. Should be indexed (.fai). "
                      "Defaults to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa"),
        ToolInput("g", Int(optional=True), prefix="-g", position=1,
                  doc="The column for gene name, or segment annotation"),
        ToolInput("H", Boolean(optional=True), prefix="-H", position=1, doc="Print this help page"),
        ToolInput("h", Boolean(optional=True), prefix="-h", position=1,
                  doc="Print a header row describing columns"),
        ToolInput("I", Int(optional=True), prefix="-I", position=1, doc="The indel size.  Default: 120bp"),
        ToolInput("i", Boolean(optional=True), prefix="-i", position=1, doc="Output splicing read counts"),
        ToolInput("k", Int(optional=True), prefix="-k", position=1,
                  doc="Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it. "
                      "For Ion or PacBio, 0 is recommended."),
        ToolInput("M", Int(optional=True), prefix="-M", position=1,
                  doc="The minimum matches for a read to be considered. If, after soft-clipping, the matched "
                      "bp is less than INT, then the read is discarded. It's meant for PCR based targeted sequencing "
                      "where there's no insert and the matching is only the primers. Default: 0, or no filtering"),
        ToolInput("m", Int(optional=True), prefix="-m", position=1,
                  doc="If set, reads with mismatches more than INT will be filtered and ignored. "
                      "Gaps are not counted as mismatches. Valid only for bowtie2/TopHat or BWA aln "
                      "followed by sampe. BWA mem is calculated as NM - Indels. "
                      "Default: 8, or reads with more than 8 mismatches will not be used."),
        ToolInput("N", String(optional=True), prefix="-N", position=1,
                  doc="The sample name to be used directly.  Will overwrite -n option"),
        ToolInput("n", String(optional=True), prefix="-n", position=1,
                  doc="The regular expression to extract sample name from BAM filenames. "
                      "Default to: /([^\/\._]+?)_[^\/]*.bam/"),
        ToolInput("O", Float(optional=True), prefix="-O", position=1,
                  doc="The reads should have at least mean MapQ to be considered a valid variant. "
                      "Default: no filtering"),
        ToolInput("o", Float(optional=True), prefix="-o", position=1,
                  doc="The Qratio of (good_quality_reads)/(bad_quality_reads+0.5). "
                      "The quality is defined by -q option.  Default: 1.5"),
        ToolInput("P", Float(optional=True), prefix="-P", position=1,
                  doc="The read position filter. If the mean variants position is less that specified, "
                      "it's considered false positive.  Default: 5"),
        ToolInput("p", Boolean(optional=True), prefix="-p", position=1,
                  doc="Do pileup regardless of the frequency"),
        ToolInput("Q", Int(optional=True), prefix="-Q", position=1,
                  doc="If set, reads with mapping quality less than INT will be filtered and ignored"),
        ToolInput("q", Int(optional=True), prefix="-q", position=1,
                  doc="The phred score for a base to be considered a good call.  "
                      "Default: 25 (for Illumina) For PGM, set it to ~15, as PGM tends to under estimate base quality."),
        ToolInput("R", String(optional=True), prefix="-R", position=1,
                  doc="The region of interest.  In the format of chr:start-end.  If end is omitted, "
                      "then a single position.  No BED is needed."),
        ToolInput("r", Int(optional=True), prefix="-r", position=1,
                  doc="The minimum # of variant reads, default 2"),
        ToolInput("S", Int(optional=True), prefix="-S", position=1,
                  doc="The column for region start, e.g. gene start"),
        ToolInput("s", Int(optional=True), prefix="-s", position=1,
                  doc="The column for segment starts in the region, e.g. exon starts"),
        ToolInput("T", Int(optional=True), prefix="-T", position=1,
                  doc="Trim bases after [INT] bases in the reads"),
        ToolInput("t", Boolean(optional=True), prefix="-t", position=1,
                  doc="Indicate to remove duplicated reads.  Only one pair with same start positions will be kept"),
        ToolInput("th", Int(optional=True), prefix="-th", position=1, doc="Threads count."),
        ToolInput("V", Int(optional=True), prefix="-V", position=1,
                  doc="The lowest frequency in the normal sample allowed for a putative somatic mutation. "
                      "Defaults to 0.05"),
        ToolInput("v", Boolean(optional=True), prefix="-v", position=1, doc="VCF format output"),
        ToolInput("VS", String(optional=True), prefix="-VS", position=1,
                  doc="[STRICT | LENIENT | SILENT] How strict to be when reading a SAM or BAM. \n"
                      "STRICT   - throw an exception if something looks wrong. "
                      "LENIENT	- Emit warnings but keep going if possible. "
                      "SILENT	- Like LENIENT, only don't emit warning messages. "
                      "Default: LENIENT"),
        ToolInput("X", Int(optional=True), prefix="-X", position=1,
                  doc="Extension of bp to look for mismatches after insersion or deletion.  "
                      "Default to 3 bp, or only calls when they're within 3 bp."),
        ToolInput("x", Int(optional=True), prefix="-x", position=1,
                  doc="The number of nucleotide to extend for each segment, default: 0"),
        ToolInput("y", Boolean(optional=True), prefix="-y", position=1, doc="<No content>"),
        ToolInput("Z", Int(optional=True), prefix="-Z", position=1,
                  doc="For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  "
                      "Default: No downsampling.  Use with caution.  "
                      "The downsampling will be random and non-reproducible."),
        ToolInput("z", Int(optional=True), prefix="-z", position=1,
                  doc="0/1  Indicate whether coordinates are zero-based, as IGV uses.  "
                      "Default: 1 for BED file or amplicon BED file. Use 0 to turn it off. "
                      "When using the -R option, it's set to 0"),
    ]

if __name__ == "__main__":
    print(SplitMultiAllele().help())
