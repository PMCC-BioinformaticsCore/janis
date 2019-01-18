from abc import ABC

from Pipeline import ToolInput, File, Boolean, String, Array, Int, Filename, ToolOutput
from Pipeline.bioinformatics.data_types.vcf import Vcf
from Pipeline.bioinformatics.tools.bcftools.bcftoolstoolbase import BcfToolsToolBase


class BcfToolsAnnotateBase(BcfToolsToolBase, ABC):

    @staticmethod
    def tool():
        return "bcftoolsAnnotate"

    @classmethod
    def bcftools_command(cls):
        return "annotate"

    def inputs(self):
        return [
            ToolInput("file", File(), position=100),
            ToolInput("outputFilename", Filename(extension=".vcf"), prefix="--output", doc='[-o] see Common Options'),
            *self.additional_args
        ]

    def outputs(self):
        return [
            ToolOutput("output", Vcf(), glob="$(inputs.outputFilename)")
        ]

    @staticmethod
    def doc():
        return BcfToolsToolBase.doc() + """
    -------------------------------------------------------------------------

    Add or remove annotations.

    Documentation: https://samtools.github.io/bcftools/bcftools.html#annotate
    """

    additional_args = [
        ToolInput("annotations", File(optional=True), prefix="--annotations",
                  doc='[-a] Bgzip-compressed and tabix-indexed file with annotations. The file can be VCF, BED, '
                      'or a tab-delimited file with mandatory columns CHROM, POS (or, alternatively, FROM and TO), '
                      'optional columns REF and ALT, and arbitrary number of annotation columns. BED files are '
                      'expected to have the ".bed" or ".bed.gz" suffix (case-insensitive), otherwise a '
                      'tab-delimited file is assumed. Note that in case of tab-delimited file, the coordinates POS, '
                      'FROM and TO are one-based and inclusive. When REF and ALT are present, only matching VCF '
                      'records will be annotated. When multiple ALT alleles are present in the annotation file '
                      '(given as comma-separated list of alleles), at least one must match one of the alleles in '
                      'the corresponding VCF record. Similarly, at least one alternate allele from a multi-allelic '
                      'VCF record must be present in the annotation file. Missing values can be added by '
                      'providing "." in place of actual value. Note that flag types, such as "INFO/FLAG", '
                      'can be annotated by including a field with the value "1" to set the flag, "0" to remove '
                      'it, or "." to keep existing flags. See also -c, --columns and -h, --header-lines.'),
        ToolInput("collapse", String(optional=True), prefix="--collapse",
                  doc='(snps|indels|both|all|some|none) Controls how to match records from the annotation file to '
                      'the target VCF. Effective only when -a is a VCF or BCF. See Common Options for more.'),
        ToolInput("columns", Array(String(), optional=True), prefix="--columns",
                  doc='[-c] Comma-separated list of columns or tags to carry over from the annotation file '
                      '(see also -a, --annotations). If the annotation file is not a VCF/BCF, list describes '
                      'the columns of the annotation file and must include CHROM, POS '
                      '(or, alternatively, FROM and TO), and optionally REF and ALT. Unused columns which '
                      'should be ignored can be indicated by "-". If the annotation file is a VCF/BCF, '
                      'only the edited columns/tags must be present and their order does not matter. '
                      'The columns ID, QUAL, FILTER, INFO and FORMAT can be edited, where INFO tags '
                      'can be written both as "INFO/TAG" or simply "TAG", and FORMAT tags can be written '
                      'as "FORMAT/TAG" or "FMT/TAG". The imported VCF annotations can be renamed as '
                      '"DST_TAG:=SRC_TAG" or "FMT/DST_TAG:=FMT/SRC_TAG". To carry over all INFO annotations, '
                      'use "INFO". To add all INFO annotations except "TAG", use "^INFO/TAG". '
                      'By default, existing values are replaced. To add annotations without overwriting '
                      'existing values (that is, to add missing tags or add values to existing tags with '
                      'missing values), use "+TAG" instead of "TAG". To append to existing values (rather '
                      'than replacing or leaving untouched), use "=TAG" (instead of "TAG" or "+TAG"). '
                      'To replace only existing values without modifying missing annotations, use "-TAG". '
                      'If the annotation file is not a VCF/BCF, all new annotations must be '
                      'defined via -h, --header-lines.'),
        ToolInput("exclude", String(optional=True), prefix="--exclude",
                  doc='[-e] exclude sites for which EXPRESSION is true. For valid expressions see EXPRESSIONS.'),
        ToolInput("headerLines", File(optional=True), prefix="--header-lines",
                  doc='[-h] Lines to append to the VCF header, see also -c, --columns and -a, --annotations.'),
        ToolInput("setId", String(optional=True), prefix="--set-id",
                  doc='[-I] assign ID on the fly. The format is the same as in the query command (see below). '
                      'By default all existing IDs are replaced. If the format string is preceded by "+", only'
                      ' missing IDs will be set. For example, one can use '
                      '# bcftools annotate --set-id +\' % CHROM\_ % POS\_ % REF\_ % FIRST_ALT\' file.vcf'),
        ToolInput("include", String(optional=True), prefix="--include",
                  doc='[-i] include only sites for which EXPRESSION is true. For valid expressions see EXPRESSIONS.'),
        ToolInput("keepSites", Boolean(optional=True), prefix="--keep-sites",
                  doc='keep sites wich do not pass -i and -e expressions instead of discarding them('),
        ToolInput("markSites", String(optional=True), prefix="--mark-sites",
                  doc='[-m] (+|-)annotate sites which are present ("+") or absent ("-") in the -a file with a '
                      'new INFO/TAG flag'),
        ToolInput("outputType", String(optional=True), prefix="--output-type",
                  doc='[-O] (b|u|z|v) see Common Options'),
        ToolInput("regions", String(optional=True), prefix="--regions",
                  doc='([-r] chr|chr:pos|chr:from-to|chr:from-[,…]) see Common Options'),
        ToolInput("regionsFile", File(optional=True), prefix="--regions-file", doc='[-R] see Common Options'),
        ToolInput("renameChrs", File(optional=True), prefix="--rename-chrs",
                  doc='rename chromosomes according to the map in file, with "old_name new_name\\n" pairs '
                      'separated by whitespaces, each on a separate line.'),
        ToolInput("samples", Array(File(), optional=True), prefix="--samples",
                  doc='[-s] subset of samples to annotate, see also Common Options'),
        ToolInput("samplesFile", File(optional=True), prefix="--samples-file",
                  doc='[-S] subset of samples to annotate. If the samples are named differently in the target '
                      'VCF and the -a, --annotations VCF, the name mapping can be given as "src_name dst_name\\n", '
                      'separated by whitespaces, each pair on a separate line.'),
        ToolInput("threads", Int(optional=True), prefix="--threads", doc='see Common Options'),
        ToolInput("remove", Array(String(), optional=True), prefix="--remove",
                  doc='[-x] List of annotations to remove. Use "FILTER" to remove all filters or "FILTER/SomeFilter" '
                      'to remove a specific filter. Similarly, "INFO" can be used to remove all INFO tags and '
                      '"FORMAT" to remove all FORMAT tags except GT. To remove all INFO tags except "FOO" and '
                      '"BAR", use "^INFO/FOO,INFO/BAR" (and similarly for FORMAT and FILTER). "INFO" can be '
                      'abbreviated to "INF" and "FORMAT" to "FMT".'),
    ]
