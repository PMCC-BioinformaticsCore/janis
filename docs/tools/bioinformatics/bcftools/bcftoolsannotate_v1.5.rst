:orphan:


BCFTools: Annotate
=====================================

Description
-------------

Tool identifier: ``bcftoolsAnnotate``

Tool path: ``janis_bioinformatics.tools.bcftools.annotate.versions import BcfToolsAnnotate_1_5``

Version: v1.5

Container: ``biocontainers/bcftools:v1.5_cv2``

Versions
*********

- `v1.9 <bcftoolsannotate_v1.9.html>`_
- v1.5 (current)

Documentation
-------------

URL
******
`https://samtools.github.io/bcftools/bcftools.html#annotate <https://samtools.github.io/bcftools/bcftools.html#annotate>`_

Tool documentation
******************
------------------------------------

Add or remove annotations.------------------------------------

Add or remove annotations.

Outputs
-------
======  =============  ===============
name    type           documentation
======  =============  ===============
out     CompressedVCF
======  =============  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======  =============  ========  ==========  ===============
name    type           prefix      position  documentation
======  =============  ========  ==========  ===============
file    CompressedVCF                   100
======  =============  ========  ==========  ===============

Optional inputs
***************

==============  =======================  ==============  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name            type                     prefix          position    documentation
==============  =======================  ==============  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
outputFilename  Optional<Filename>       --output                    [-o] see Common Options
annotations     Optional<File>           --annotations               [-a] Bgzip-compressed and tabix-indexed file with annotations. The file can be VCF, BED, or a tab-delimited file with mandatory columns CHROM, POS (or, alternatively, FROM and TO), optional columns REF and ALT, and arbitrary number of annotation columns. BED files are expected to have the ".bed" or ".bed.gz" suffix (case-insensitive), otherwise a tab-delimited file is assumed. Note that in case of tab-delimited file, the coordinates POS, FROM and TO are one-based and inclusive. When REF and ALT are present, only matching VCF records will be annotated. When multiple ALT alleles are present in the annotation file (given as comma-separated list of alleles), at least one must match one of the alleles in the corresponding VCF record. Similarly, at least one alternate allele from a multi-allelic VCF record must be present in the annotation file. Missing values can be added by providing "." in place of actual value. Note that flag types, such as "INFO/FLAG", can be annotated by including a field with the value "1" to set the flag, "0" to remove it, or "." to keep existing flags. See also -c, --columns and -h, --header-lines.
collapse        Optional<String>         --collapse                  (snps|indels|both|all|some|none) Controls how to match records from the annotation file to the target VCF. Effective only when -a is a VCF or BCF. See Common Options for more.
columns         Optional<Array<String>>  --columns                   [-c] Comma-separated list of columns or tags to carry over from the annotation file (see also -a, --annotations). If the annotation file is not a VCF/BCF, list describes the columns of the annotation file and must include CHROM, POS (or, alternatively, FROM and TO), and optionally REF and ALT. Unused columns which should be ignored can be indicated by "-". If the annotation file is a VCF/BCF, only the edited columns/tags must be present and their order does not matter. The columns ID, QUAL, FILTER, INFO and FORMAT can be edited, where INFO tags can be written both as "INFO/TAG" or simply "TAG", and FORMAT tags can be written as "FORMAT/TAG" or "FMT/TAG". The imported VCF annotations can be renamed as "DST_TAG:=SRC_TAG" or "FMT/DST_TAG:=FMT/SRC_TAG". To carry over all INFO annotations, use "INFO". To add all INFO annotations except "TAG", use "^INFO/TAG". By default, existing values are replaced. To add annotations without overwriting existing values (that is, to add missing tags or add values to existing tags with missing values), use "+TAG" instead of "TAG". To append to existing values (rather than replacing or leaving untouched), use "=TAG" (instead of "TAG" or "+TAG"). To replace only existing values without modifying missing annotations, use "-TAG". If the annotation file is not a VCF/BCF, all new annotations must be defined via -h, --header-lines.
exclude         Optional<String>         --exclude                   [-e] exclude sites for which EXPRESSION is true. For valid expressions see EXPRESSIONS.
headerLines     Optional<File>           --header-lines              [-h] Lines to append to the VCF header, see also -c, --columns and -a, --annotations.
setId           Optional<String>         --set-id                    [-I] assign ID on the fly. The format is the same as in the query command (see below). By default all existing IDs are replaced. If the format string is preceded by "+", only missing IDs will be set. For example, one can use # bcftools annotate --set-id +' % CHROM\_ % POS\_ % REF\_ % FIRST_ALT' file.vcf
include         Optional<String>         --include                   [-i] include only sites for which EXPRESSION is true. For valid expressions see EXPRESSIONS.
keepSites       Optional<Boolean>        --keep-sites                keep sites wich do not pass -i and -e expressions instead of discarding them(
markSites       Optional<String>         --mark-sites                [-m] (+|-)annotate sites which are present ("+") or absent ("-") in the -a file with a new INFO/TAG flag
outputType      Optional<String>         --output-type               [-O] (b|u|z|v) see Common Options
regions         Optional<String>         --regions                   ([-r] chr|chr:pos|chr:from-to|chr:from-[,…]) see Common Options
regionsFile     Optional<File>           --regions-file              [-R] see Common Options
renameChrs      Optional<File>           --rename-chrs               rename chromosomes according to the map in file, with "old_name new_name\n" pairs separated by whitespaces, each on a separate line.
samples         Optional<Array<File>>    --samples                   [-s] subset of samples to annotate, see also Common Options
samplesFile     Optional<File>           --samples-file              [-S] subset of samples to annotate. If the samples are named differently in the target VCF and the -a, --annotations VCF, the name mapping can be given as "src_name dst_name\n", separated by whitespaces, each pair on a separate line.
threads         Optional<Integer>        --threads                   see Common Options
remove          Optional<Array<String>>  --remove                    [-x] List of annotations to remove. Use "FILTER" to remove all filters or "FILTER/SomeFilter" to remove a specific filter. Similarly, "INFO" can be used to remove all INFO tags and "FORMAT" to remove all FORMAT tags except GT. To remove all INFO tags except "FOO" and "BAR", use "^INFO/FOO,INFO/BAR" (and similarly for FORMAT and FILTER). "INFO" can be abbreviated to "INF" and "FORMAT" to "FMT".
==============  =======================  ==============  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: **Unknown**


*BCFTools: Annotate was last updated on 2019-01-24*.
*This page was automatically generated on 2019-09-26*.
