:orphan:


BCFTools: Normalize
==================================

Description
-------------

Tool identifier: ``bcftoolsNorm``

Tool path: ``janis_bioinformatics.tools.bcftools.norm.norm_1_5 import BcfToolsNorm_1_5``

Version: v1.5

Docker: ``biocontainers/bcftools:v1.5_cv2``

Versions
*********

- `v1.9 <bcftoolsnorm_v1.9.html>`_
- v1.5 (current)

Documentation
-------------

URL
******
*No URL to the documentation was provided*

Tool documentation
******************
*No documentation was provided: `contribute one <https://github.com/illusional>`_*

Outputs
-------
======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======  ======  ========  ==========  ===============
name    type    prefix      position  documentation
======  ======  ========  ==========  ===============
vcf     VCF                       10
======  ======  ========  ==========  ===============

Optional inputs
***************

=====================  =======================  ============  ==========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                   type                     prefix        position    documentation
=====================  =======================  ============  ==========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
outputFilename         Optional<Filename>       -o                        --output: When output consists of a single stream, write it to FILE rather than to standard output, where it is written by default.
checkRef               Optional<String>         -c                        --check-ref e|w|x|s: what to do when incorrect or missing REF allele is encountered: exit (e), warn (w), exclude (x), or set/fix (s) bad sites. The w option can be combined with x and s. Note that s can swap alleles and will update genotypes (GT) and AC counts, but will not attempt to fix PL or other fields. Also note, and this cannot be stressed enough, that s will NOT fix strand issues in your VCF, do NOT use it for that purpose!!! (Instead see http://samtools.github.io/bcftools/howtos/plugin.af-dist.html and http://samtools.github.io/bcftools/howtos/plugin.fixref.html.)
removeDups             Optional<String>         -d                        --rm-dup: snps|indels|both|all|none. If a record is present multiple times, output only the first instance, see --collapse in Common Options.
removeDupsAcrossFiles  Optional<Boolean>        -D                        --remove-duplicates: If a record is present in multiple files, output only the first instance. Alias for -d none, deprecated.
reference              Optional<FastaWithDict>  -f                        --fasta-ref: reference sequence. Supplying this option will turn on left-alignment and normalization, however, see also the --do-not-normalize option below.
multiallelics          Optional<String>         -m                        --multiallelics -|+[snps|indels|both|any]: split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+). An optional type string can follow which controls variant types which should be split or merged together: If only SNP records should be split or merged, specify snps; if both SNPs and indels should be merged separately into two records, specify both; if SNPs and indels should be merged into a single record, specify any.
noVersion              Optional<Boolean>        --no-version              Do not append version and command line information to the output VCF header.
noNormalize            Optional<Boolean>        -N                        --do-not-normalize: the -c s option can be used to fix or set the REF allele from the reference -f. The -N option will not turn on indel normalisation as the -f option normally implies
outputType             Optional<String>         -O                        --output-type b|u|z|v: Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion.
regions                Optional<String>         -r                        --regions chr|chr:pos|chr:from-to|chr:from-[,…]: Comma-separated list of regions, see also -R, --regions-file. Note that -r cannot be used in combination with -R.
regionsFile            Optional<File>           -R                        --regions-file: Regions can be specified either on command line or in a VCF, BED, or tab-delimited file (the default). The columns of the tab-delimited file are: CHROM, POS, and, optionally, POS_TO, where positions are 1-based and inclusive. The columns of the tab-delimited BED file are also CHROM, POS and POS_TO (trailing columns are ignored), but coordinates are 0-based, half-open. To indicate that a file be treated as BED rather than the 1-based tab-delimited file, the file must have the '.bed' or '.bed.gz' suffix (case-insensitive). Uncompressed files are stored in memory, while bgzip-compressed and tabix-indexed region files are streamed. Note that sequence names must match exactly, 'chr20' is not the same as '20'. Also note that chromosome ordering in FILE will be respected, the VCF will be processed in the order in which chromosomes first appear in FILE. However, within chromosomes, the VCF will always be processed in ascending genomic coordinate order no matter what order they appear in FILE. Note that overlapping regions in FILE can result in duplicated out of order positions in the output. This option requires indexed VCF/BCF files. Note that -R cannot be used in combination with -r.
strictFilter           Optional<Boolean>        -s                        --strict-filter: when merging (-m+), merged site is PASS only if all sites being merged PASS
targets                Optional<Array<File>>    -t                        --targets: [^]chr|chr:pos|chr:from-to|chr:from-[,…]: Similar as -r, --regions, but the next position is accessed by streaming the whole VCF/BCF rather than using the tbi/csi index. Both -r and -t options can be applied simultaneously: -r uses the index to jump to a region and -t discards positions which are not in the targets. Unlike -r, targets can be prefixed with '^' to request logical complement. For example, '^X,Y,MT' indicates that sequences X, Y and MT should be skipped. Yet another difference between the two is that -r checks both start and end positions of indels, whereas -t checks start positions only. Note that -t cannot be used in combination with -T.
targetsFile            Optional<File>           -T                        --targets-file: Same -t, --targets, but reads regions from a file. Note that -T cannot be used in combination with -t. With the call -C alleles command, third column of the targets file must be comma-separated list of alleles, starting with the reference allele. Note that the file must be compressed and index. Such a file can be easily created from a VCF using: `bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' file.vcf | bgzip -c > als.tsv.gz && tabix -s1 -b2 -e2 als.tsv.gz`
threads                Optional<Integer>        --threads                 Number of output compression threads to use in addition to main thread. Only used when --output-type is b or z. Default: 0.
siteWin                Optional<Integer>        -w                        --site-win: maximum distance between two records to consider when locally sorting variants which changed position during the realignment
=====================  =======================  ============  ==========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: **Unknown**


*BCFTools: Normalize was last updated on **Unknown***.
*This page was automatically generated on 2019-07-26*.
