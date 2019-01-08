from abc import ABC

from Pipeline import ToolInput, String, Boolean, File, Filename, Array, Int, ToolOutput
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.vcf import TabixIdx, Vcf, VcfIdx
from Pipeline.bioinformatics.tools.bcftools.bcftoolstoolbase import BcfToolsToolBase


class BcfToolsNormBase(BcfToolsToolBase):

    @staticmethod
    def tool():
        return "bcftools-norm"

    @classmethod
    def bcftools_command(cls):
        return "norm"

    def inputs(self):
        return [
            ToolInput("input", Vcf(), position=10),
            ToolInput("outputFilename", Filename(extension=".vcf"), prefix="-o",
                      doc="--output: When output consists of a single stream, "
                          "write it to FILE rather than to standard output, where it is written by default."),
            *self.additional_args
        ]

    def outputs(self):
        return [
            ToolOutput("output", Vcf(), glob="$(inputs.outputFilename)")
        ]

    @staticmethod
    def doc():
        return BcfToolsToolBase.doc() + """\n
    -------------------------------------------------------------------------
    
    Left-align and normalize indels, check if REF alleles match the reference, 
    split multiallelic sites into multiple rows; recover multiallelics from 
    multiple rows. Left-alignment and normalization will only be applied if 
    the --fasta-ref option is supplied.
    
    Documentation: https://samtools.github.io/bcftools/bcftools.html#norm
    """.strip()

    additional_args = [

        ToolInput("checkRef", String(optional=True), prefix="-c",
                  doc="--check-ref e|w|x|s: what to do when incorrect or missing REF allele is encountered: "
                      "exit (e), warn (w), exclude (x), or set/fix (s) bad sites. The w option can be combined with "
                      "x and s. Note that s can swap alleles and will update genotypes (GT) and AC counts, but will "
                      "not attempt to fix PL or other fields. Also note, and this cannot be stressed enough, that s "
                      "will NOT fix strand issues in your VCF, do NOT use it for that purpose!!! "
                      "(Instead see http://samtools.github.io/bcftools/howtos/plugin.af-dist.html "
                      "and http://samtools.github.io/bcftools/howtos/plugin.fixref.html.)"),
        ToolInput("removeDups", String(optional=True), prefix="-d",
                  doc="--rm-dup: snps|indels|both|all|none. If a record is present multiple times, "
                      "output only the first instance, see --collapse in Common Options."),
        ToolInput("removeDupsAcrossFiles", Boolean(optional=True), prefix="-D",
                  doc="--remove-duplicates: If a record is present in multiple files, "
                      "output only the first instance. Alias for -d none, deprecated."),
        ToolInput("reference", FastaWithDict(optional=True), prefix="-f",
                  doc="--fasta-ref: reference sequence. Supplying this option will turn on left-alignment and "
                      "normalization, however, see also the --do-not-normalize option below."),
        ToolInput("multiallelics", String(optional=True), prefix="-m",
                  doc="--multiallelics -|+[snps|indels|both|any]: split multiallelic sites into "
                      "biallelic records (-) or join biallelic sites into multiallelic records (+). "
                      "An optional type string can follow which controls variant types which should "
                      "be split or merged together: If only SNP records should be split or merged, "
                      "specify snps; if both SNPs and indels should be merged separately into two "
                      "records, specify both; if SNPs and indels should be merged into a single record, specify any."),
        ToolInput("noVersion", Boolean(optional=True), prefix="--no-version",
                  doc="Do not append version and command line information to the output VCF header."),
        ToolInput("noNormalize", Boolean(optional=True), prefix="-N",
                  doc="--do-not-normalize: the -c s option can be used to fix or set the REF allele from the reference"
                      " -f. The -N option will not turn on indel normalisation as the -f option normally implies"),
        ToolInput("outputType", String(optional=True), prefix="-O",
                  doc="--output-type b|u|z|v: Output compressed BCF (b), uncompressed BCF (u), "
                      "compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping "
                      "between bcftools subcommands to speed up performance by removing "
                      "unnecessary compression/decompression and VCF←→BCF conversion."),
        ToolInput("regions", String(optional=True), prefix="-r",
                  doc="--regions chr|chr:pos|chr:from-to|chr:from-[,…]: Comma-separated list of regions, "
                      "see also -R, --regions-file. Note that -r cannot be used in combination with -R."),
        ToolInput("regionsFile", File(optional=True), prefix="-R",
                  doc="--regions-file: Regions can be specified either on command line or in a VCF, BED, or "
                      "tab-delimited file (the default). The columns of the tab-delimited file are: CHROM, POS, "
                      "and, optionally, POS_TO, where positions are 1-based and inclusive. The columns of the "
                      "tab-delimited BED file are also CHROM, POS and POS_TO (trailing columns are ignored), "
                      "but coordinates are 0-based, half-open. To indicate that a file be treated as BED rather "
                      "than the 1-based tab-delimited file, the file must have the '.bed' or '.bed.gz' suffix "
                      "(case-insensitive). Uncompressed files are stored in memory, while bgzip-compressed and "
                      "tabix-indexed region files are streamed. Note that sequence names must match exactly, 'chr20'"
                      " is not the same as '20'. Also note that chromosome ordering in FILE will be respected, "
                      "the VCF will be processed in the order in which chromosomes first appear in FILE. "
                      "However, within chromosomes, the VCF will always be processed in ascending genomic coordinate "
                      "order no matter what order they appear in FILE. Note that overlapping regions in FILE can "
                      "result in duplicated out of order positions in the output. This option requires indexed "
                      "VCF/BCF files. Note that -R cannot be used in combination with -r."),
        ToolInput("strictFilter", Boolean(optional=True), prefix="-s",
                  doc="--strict-filter: when merging (-m+), merged site is PASS only if all sites being merged PASS"),
        ToolInput("targets", Array(File(), optional=True), prefix="-t",
                  doc="--targets: [^]chr|chr:pos|chr:from-to|chr:from-[,…]: Similar as -r, --regions, but the next "
                      "position is accessed by streaming the whole VCF/BCF rather than using the tbi/csi index. "
                      "Both -r and -t options can be applied simultaneously: -r uses the index to jump to a region "
                      "and -t discards positions which are not in the targets. Unlike -r, targets can be prefixed "
                      "with '^' to request logical complement. For example, '^X,Y,MT' indicates that sequences "
                      "X, Y and MT should be skipped. Yet another difference between the two is that -r checks "
                      "both start and end positions of indels, whereas -t checks start positions only. "
                      "Note that -t cannot be used in combination with -T. "),
        ToolInput("targetsFile", File(optional=True), prefix="-T",
                  doc="--targets-file: Same -t, --targets, but reads regions from a file. "
                      "Note that -T cannot be used in combination with -t. With the call -C alleles command, "
                      "third column of the targets file must be comma-separated list of alleles, starting with "
                      "the reference allele. Note that the file must be compressed and index. "
                      "Such a file can be easily created from a VCF using: "
                      "`bcftools query -f'%CHROM\\t%POS\\t%REF,%ALT\\n' file.vcf | bgzip -c > als.tsv.gz "
                      "&& tabix -s1 -b2 -e2 als.tsv.gz`"),
        ToolInput("threads", Int(optional=True), prefix="--threads",
                  doc="Number of output compression threads to use in addition to main thread. "
                      "Only used when --output-type is b or z. Default: 0."),
        ToolInput("siteWin", Int(optional=True), prefix="-w",
                  doc="--site-win: maximum distance between two records to consider when locally "
                      "sorting variants which changed position during the realignment")
    ]
