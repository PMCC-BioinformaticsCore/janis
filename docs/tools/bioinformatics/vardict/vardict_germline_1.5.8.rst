:orphan:

VarDict (Germline)
=====================================

0 contributors · 5 versions

:ID: ``vardict_germline``
:Python: ``janis_bioinformatics.tools.vardict.vardictgermline import VarDictGermline_1_5_8``
:Versions: 1.7.0, 1.6.0, 1.5.8, 1.5.7, 1.5.6
:Container: michaelfranklin/vardict:1.5.8
:Authors: 
:Citations: None
:Created: None
:Updated: None
:Required inputs:
   - ``intervals: bed``

   - ``bam: BamPair``

   - ``reference: FastaFai``

   - ``sampleName: String``

   - ``var2vcfSampleName: String``

   - ``var2vcfAlleleFreqThreshold: Float``
:Outputs: 
   - ``out: VCF``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

------

Additional configuration (inputs)
---------------------------------

==========================  ==================  ==================================================================================================================================================================================================================================================================================
name                        type                documentation
==========================  ==================  ==================================================================================================================================================================================================================================================================================
intervals                   bed
bam                         BamPair             The indexed BAM file
reference                   FastaFai            The reference fasta. Should be indexed (.fai). Defaults to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa
sampleName                  String              The sample name to be used directly.  Will overwrite -n option
var2vcfSampleName           String
var2vcfAlleleFreqThreshold  Float
outputFilename              Optional<Filename>
indels3prime                Optional<Boolean>   Indicate to move indels to 3-prime if alternative alignment can be achieved.
amplicon                    Optional<Float>     Indicate it's amplicon based calling.  Reads that don't map to the amplicon will be skipped.  A read pair is considered belonging  to the amplicon if the edges are less than int bp to the amplicon, and overlap fraction is at least float.  Default: 10:0.95
minReads                    Optional<Integer>   The minimum # of reads to determine strand bias, default 2
chromNamesAreNumbers        Optional<Boolean>   Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
chromColumn                 Optional<Integer>   The column for chromosome
debug                       Optional<Boolean>   Debug mode.  Will print some error messages and append full genotype at the end.
splitDelimeter              Optional<String>    The delimiter for split region_info, default to tab "	"
geneEndCol                  Optional<Integer>   The column for region end, e.g. gene end
segEndCol                   Optional<Integer>   The column for segment ends in the region, e.g. exon ends
filter                      Optional<String>    The hexical to filter reads using samtools. Default: 0x500 (filter 2nd alignments and duplicates). Use -F 0 to turn it off.
alleleFreqThreshold         Optional<Float>     The threshold for allele frequency, default: 0.05 or 5%
geneNameCol                 Optional<Integer>   The column for gene name, or segment annotation
printHeaderRow              Optional<Boolean>   Print a header row describing columns
indelSize                   Optional<Integer>   The indel size.  Default: 120bp
outputSplice                Optional<Boolean>   Output splicing read counts
performLocalRealignment     Optional<Integer>   Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it. For Ion or PacBio, 0 is recommended.
minMatches                  Optional<Integer>   The minimum matches for a read to be considered. If, after soft-clipping, the matched bp is less than INT, then the read is discarded. It's meant for PCR based targeted sequencing where there's no insert and the matching is only the primers. Default: 0, or no filtering
maxMismatches               Optional<Integer>   If set, reads with mismatches more than INT will be filtered and ignored. Gaps are not counted as mismatches. Valid only for bowtie2/TopHat or BWA aln followed by sampe. BWA mem is calculated as NM - Indels. Default: 8, or reads with more than 8 mismatches will not be used.
regexSampleName             Optional<String>    The regular expression to extract sample name from BAM filenames. Default to: /([^\/\._]+?)_[^\/]*.bam/
mapq                        Optional<String>    The reads should have at least mean MapQ to be considered a valid variant. Default: no filtering
qratio                      Optional<Float>     The Qratio of (good_quality_reads)/(bad_quality_reads+0.5). The quality is defined by -q option.  Default: 1.5
readPosition                Optional<Float>     The read position filter. If the mean variants position is less that specified, it's considered false positive.  Default: 5
pileup                      Optional<Boolean>   Do pileup regardless of the frequency
minMappingQual              Optional<Integer>   If set, reads with mapping quality less than INT will be filtered and ignored
phredScore                  Optional<Integer>   The phred score for a base to be considered a good call.  Default: 25 (for Illumina) For PGM, set it to ~15, as PGM tends to under estimate base quality.
region                      Optional<String>    The region of interest.  In the format of chr:start-end.  If end is omitted, then a single position.  No BED is needed.
minVariantReads             Optional<Integer>   The minimum # of variant reads, default 2
regStartCol                 Optional<Integer>   The column for region start, e.g. gene start
segStartCol                 Optional<Integer>   The column for segment starts in the region, e.g. exon starts
minReadsBeforeTrim          Optional<Integer>   Trim bases after [INT] bases in the reads
removeDuplicateReads        Optional<Boolean>   Indicate to remove duplicated reads.  Only one pair with same start positions will be kept
threads                     Optional<Integer>   Threads count.
freq                        Optional<Integer>   The lowest frequency in the normal sample allowed for a putative somatic mutation. Defaults to 0.05
vcfFormat                   Optional<Boolean>   VCF format output
vs                          Optional<String>    [STRICT | LENIENT | SILENT] How strict to be when reading a SAM or BAM: STRICT   - throw an exception if something looks wrong. LENIENT	- Emit warnings but keep going if possible. SILENT	- Like LENIENT, only don't emit warning messages. Default: LENIENT
bp                          Optional<Integer>   Extension of bp to look for mismatches after insersion or deletion.  Default to 3 bp, or only calls when they're within 3 bp.
extensionNucleotide         Optional<Integer>   The number of nucleotide to extend for each segment, default: 0
yy                          Optional<Boolean>   <No content>
downsamplingFraction        Optional<Integer>   For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  Default: No downsampling.  Use with caution.  The downsampling will be random and non-reproducible.
zeroBasedCoords             Optional<Integer>   0/1  Indicate whether coordinates are zero-based, as IGV uses.  Default: 1 for BED file or amplicon BED file. Use 0 to turn it off. When using the -R option, it's set to 0
==========================  ==================  ==================================================================================================================================================================================================================================================================================

