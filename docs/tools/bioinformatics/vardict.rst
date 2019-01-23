
vardict
=======
*bioinformatics*

Documentation
-------------

Docker
******
``ubuntu:latest``

URL
******
`https://github.com/AstraZeneca-NGS/VarDict <https://github.com/AstraZeneca-NGS/VarDict>`_

Docstring
*********

    VarDict
    
    VarDict is an ultra sensitive variant caller for both single and paired sample variant 
    calling from BAM files. VarDict implements several novel features such as amplicon bias 
    aware variant calling from targeted sequencing experiments, rescue of long indels by 
    realigning bwa soft clipped reads and better scalability than many Java based variant callers.

    Due to the philosophy of VarDict in calling "everything", several downstream strategies have 
    been developed to filter variants to for example the most likely cancer driving events. 
    These strategies are based on evidence in different databases and/or quality metrics. 
    http://bcb.io/2016/04/04/vardict-filtering/ provides an overview of how to develop further 
    filters for VarDict. The script at https://github.com/AstraZeneca-NGS/VarDict/blob/master/vcf2txt.pl 
    can be used to put the variants into a context by including information from dbSNP, Cosmic and ClinVar. 
    We are open to suggestions from the community on how to best narrow down to the variants of most interest.
    
    A Java based drop-in replacement for vardict.pl is being developed at 
    https://github.com/AstraZeneca-NGS/VarDictJava. The Java implementation is approximately 
    10 times faster than the original Perl implementation and does not depend on samtools
    
    To enable amplicon aware variant calling (single sample mode only; not supported in paired 
    variant calling), please make sure the bed file has 8 columns with the 7th and 8th columns 
    containing the insert interval (therefore subset of the 2nd and 3rd column interval). 
    
    Requirements

        - Perl (uses /usr/bin/env perl)
        - R (uses /usr/bin/env R)
        - samtools (must be in path, not required if using the Java implementation in place of vardict.pl)
    

Outputs
-------
======  ======  ===============
name    type    documentation
======  ======  ===============
output  VCF
======  ======  ===============

Inputs
------
==========================  ==================  ========  ==========  ==================================================================================================================================================================================================================================================================================
name                        type                prefix      position  documentation
==========================  ==================  ========  ==========  ==================================================================================================================================================================================================================================================================================
input                       bed                                    2
outputFilename              Optional<Filename>  >                  6
indexedBam                  BamPair             -b                 1  The indexed BAM file
referenceFasta              FastaFai            -G                 1  The reference fasta. Should be indexed (.fai). Defaults to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa
indels3prime                Optional<Boolean>   -3                 1  Indicate to move indels to 3-prime if alternative alignment can be achieved.
amplicon                    Optional<Float>     -a                 1  Indicate it's amplicon based calling.  Reads that don't map to the amplicon will be skipped.  A read pair is considered belonging  to the amplicon if the edges are less than int bp to the amplicon, and overlap fraction is at least float.  Default: 10:0.95
minReads                    Optional<Integer>   -B                 1  The minimum # of reads to determine strand bias, default 2
chromNamesAreNumbers        Optional<Boolean>   -C                 1  Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
chromColumn                 Optional<Integer>   -c                 1  The column for chromosome
debug                       Optional<Boolean>   -D                 1  Debug mode.  Will print some error messages and append full genotype at the end.
splitDelimeter              Optional<String>    -d                 1  The delimiter for split region_info, default to tab "	"
geneEndCol                  Optional<Integer>   -E                 1  The column for region end, e.g. gene end
segEndCol                   Optional<Integer>   -e                 1  The column for segment ends in the region, e.g. exon ends
filter                      Optional<String>    -F                 1  The hexical to filter reads using samtools. Default: 0x500 (filter 2nd alignments and duplicates). Use -F 0 to turn it off.
alleleFreqThreshold         Optional<Float>     -f                 1  The threshold for allele frequency, default: 0.05 or 5%
geneNameCol                 Optional<Integer>   -g                 1  The column for gene name, or segment annotation
printHeaderRow              Optional<Boolean>   -h                 1  Print a header row describing columns
indelSize                   Optional<Integer>   -I                 1  The indel size.  Default: 120bp
outputSplice                Optional<Boolean>   -i                 1  Output splicing read counts
performLocalRealignment     Optional<Integer>   -k                 1  Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it. For Ion or PacBio, 0 is recommended.
minMatches                  Optional<Integer>   -M                 1  The minimum matches for a read to be considered. If, after soft-clipping, the matched bp is less than INT, then the read is discarded. It's meant for PCR based targeted sequencing where there's no insert and the matching is only the primers. Default: 0, or no filtering
maxMismatches               Optional<Integer>   -m                 1  If set, reads with mismatches more than INT will be filtered and ignored. Gaps are not counted as mismatches. Valid only for bowtie2/TopHat or BWA aln followed by sampe. BWA mem is calculated as NM - Indels. Default: 8, or reads with more than 8 mismatches will not be used.
sampleName                  String              -N                 1  The sample name to be used directly.  Will overwrite -n option
regexSampleName             Optional<String>    -n                 1  The regular expression to extract sample name from BAM filenames. Default to: /([^\/\._]+?)_[^\/]*.bam/
mapq                        Optional<String>    -O                 1  The reads should have at least mean MapQ to be considered a valid variant. Default: no filtering
qratio                      Optional<Float>     -o                 1  The Qratio of (good_quality_reads)/(bad_quality_reads+0.5). The quality is defined by -q option.  Default: 1.5
readPosition                Optional<Float>     -P                 1  The read position filter. If the mean variants position is less that specified, it's considered false positive.  Default: 5
pileup                      Optional<Boolean>   -p                 1  Do pileup regardless of the frequency
minMappingQual              Optional<Integer>   -Q                 1  If set, reads with mapping quality less than INT will be filtered and ignored
phredScore                  Optional<Integer>   -q                 1  The phred score for a base to be considered a good call.  Default: 25 (for Illumina) For PGM, set it to ~15, as PGM tends to under estimate base quality.
region                      Optional<String>    -R                 1  The region of interest.  In the format of chr:start-end.  If end is omitted, then a single position.  No BED is needed.
minVariantReads             Optional<Integer>   -r                 1  The minimum # of variant reads, default 2
regStartCol                 Optional<Integer>   -S                 1  The column for region start, e.g. gene start
segStartCol                 Optional<Integer>   -s                 1  The column for segment starts in the region, e.g. exon starts
minReadsBeforeTrim          Optional<Integer>   -T                 1  Trim bases after [INT] bases in the reads
removeDuplicateReads        Optional<Boolean>   -t                 1  Indicate to remove duplicated reads.  Only one pair with same start positions will be kept
threads                     Optional<Integer>   -th                1  Threads count.
freq                        Optional<Integer>   -V                 1  The lowest frequency in the normal sample allowed for a putative somatic mutation. Defaults to 0.05
vcfFormat                   Optional<Boolean>   -v                 1  VCF format output
vs                          Optional<String>    -VS                1  [STRICT | LENIENT | SILENT] How strict to be when reading a SAM or BAM: STRICT   - throw an exception if something looks wrong. LENIENT	- Emit warnings but keep going if possible. SILENT	- Like LENIENT, only don't emit warning messages. Default: LENIENT
bp                          Optional<Integer>   -X                 1  Extension of bp to look for mismatches after insersion or deletion.  Default to 3 bp, or only calls when they're within 3 bp.
extensionNucleotide         Optional<Integer>   -x                 1  The number of nucleotide to extend for each segment, default: 0
y                           Optional<Boolean>   -y                 1  <No content>
downsamplingFraction        Optional<Integer>   -Z                 1  For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  Default: No downsampling.  Use with caution.  The downsampling will be random and non-reproducible.
zeroBasedCoords             Optional<Integer>   -z                 1  0/1  Indicate whether coordinates are zero-based, as IGV uses.  Default: 1 for BED file or amplicon BED file. Use 0 to turn it off. When using the -R option, it's set to 0
var2vcfSampleName           String              -N                 5
var2vcfAlleleFreqThreshold  Float               -f                 5
==========================  ==================  ========  ==========  ==================================================================================================================================================================================================================================================================================


*This page was automatically generated*
