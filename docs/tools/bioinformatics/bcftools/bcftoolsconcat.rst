
.. include:: bcftoolsconcat_v1.9

BCFTools: Concat
=================================

Description
-------------

Tool identifier: ``bcftoolsConcat``

Tool path: ``janis_bioinformatics.tools.bcftools.concat.versions import BcfToolsConcat_1_9``

Version: v1.9

Container: ``michaelfranklin/bcftools:1.9``



Documentation
-------------

URL
******
`https://samtools.github.io/bcftools/bcftools.html#concat <https://samtools.github.io/bcftools/bcftools.html#concat>`_

Tool documentation
******************

Concatenate or combine VCF/BCF files. All source files must have the same sample
columns appearing in the same order. The program can be used, for example, to
concatenate chromosome VCFs into one VCF, or combine a SNP VCF and an indel
VCF into one. The input files must be sorted by chr and position. The files
must be given in the correct order to produce sorted VCF on output unless
the -a, --allow-overlaps option is specified. With the --naive option, the files
are concatenated without being recompressed, which is very fast but dangerous
if the BCF headers differ.


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

======  ====================  ========  ==========  ===============
name    type                  prefix      position  documentation
======  ====================  ========  ==========  ===============
vcf     Array<CompressedVCF>                    15
======  ====================  ========  ==========  ===============

Optional inputs
***************

==============  ==================  ============  ==========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name            type                prefix        position    documentation
==============  ==================  ============  ==========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
outputFilename  Optional<Filename>  -o                        --output: When output consists of a single stream, write it to FILE rather than to standard output, where it is written by default.
allowOverLaps   Optional<Boolean>   -a                        First coordinate of the next file can precede last record of the current file.
compactPS       Optional<Boolean>   -c                        Do not output PS tag at each site, only at the start of a new phase set block.
rmDups          Optional<String>    -d                        Output duplicate records present in multiple files only once: <snps|indels|both|all|none>
rmDupsNone      Optional<Boolean>   -d                        Alias for -d none
fileList        Optional<File>      -f                        Read the list of files from a file.
ligate          Optional<Boolean>   -l                        Ligate phased VCFs by matching phase at overlapping haplotypes
noVersion       Optional<Boolean>   --no-version              Do not append version and command line information to the output VCF header.
naive           Optional<Boolean>   -n                        Concatenate files without recompression (dangerous, use with caution)
outputType      Optional<String>    -O                        --output-type b|u|z|v: Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion.
minPG           Optional<Integer>   -q                        Break phase set if phasing quality is lower than <int> [30]
regions         Optional<String>    -r                        --regions chr|chr:pos|chr:from-to|chr:from-[,…]: Comma-separated list of regions, see also -R, --regions-file. Note that -r cannot be used in combination with -R.
regionsFile     Optional<File>      -R                        --regions-file: Regions can be specified either on command line or in a VCF, BED, or tab-delimited file (the default). The columns of the tab-delimited file are: CHROM, POS, and, optionally, POS_TO, where positions are 1-based and inclusive. The columns of the tab-delimited BED file are also CHROM, POS and POS_TO (trailing columns are ignored), but coordinates are 0-based, half-open. To indicate that a file be treated as BED rather than the 1-based tab-delimited file, the file must have the '.bed' or '.bed.gz' suffix (case-insensitive). Uncompressed files are stored in memory, while bgzip-compressed and tabix-indexed region files are streamed. Note that sequence names must match exactly, 'chr20' is not the same as '20'. Also note that chromosome ordering in FILE will be respected, the VCF will be processed in the order in which chromosomes first appear in FILE. However, within chromosomes, the VCF will always be processed in ascending genomic coordinate order no matter what order they appear in FILE. Note that overlapping regions in FILE can result in duplicated out of order positions in the output. This option requires indexed VCF/BCF files. Note that -R cannot be used in combination with -r.
threads         Optional<Integer>   --threads                 Number of output compression threads to use in addition to main thread. Only used when --output-type is b or z. Default: 0.
==============  ==================  ============  ==========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: **Unknown**


*BCFTools: Concat was last updated on 2019-09-09*.
*This page was automatically generated on 2019-09-26*.
