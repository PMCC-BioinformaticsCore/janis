:orphan:

BCFTools: View
=============================

0 contributors · 2 versions

:ID: ``bcftoolsview``
:Python: ``janis_bioinformatics.tools.bcftools.view.versions import BcfToolsView_1_5``
:Versions: v1.9, v1.5
:Container: biocontainers/bcftools:v1.5_cv2
:Authors: 
:Citations: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9
:DOI: http://www.ncbi.nlm.nih.gov/pubmed/19505943
:Created: None
:Updated: 2019-01-24
:Required inputs:
   - ``file: CompressedVCF``
:Outputs: 
   - ``out: stdout<CompressedVCF>``

Documentation
-------------

URL: `https://samtools.github.io/bcftools/bcftools.html#view <https://samtools.github.io/bcftools/bcftools.html#view>`_

________________________________
 
        View, subset and filter VCF or BCF files by position and filtering expression
        Convert between VCF and BCF. Former bcftools subset.

------

Additional configuration (inputs)
---------------------------------

================  =======================  ==============================================================================================================================================================================
name              type                     documentation
================  =======================  ==============================================================================================================================================================================
file              CompressedVCF
dropGenotypes     Optional<Boolean>        (-G) drop individual genotype information (after subsetting if -s option set)
headerOnly        Optional<Boolean>        (-h) print the header only
noHeader          Optional<Boolean>        (-H) suppress the header in VCF output
compressionLevel  Optional<Integer>        (-l) compression level: 0 uncompressed, 1 best speed, 9 best compression [-1]
noVersion         Optional<Boolean>        do not append version and command line to the header
regions           Optional<String>         (-r) restrict to comma-separated list of regions
regionsFile       Optional<File>           (-R) restrict to regions listed in a file
targets           Optional<String>         (-t) similar to -r but streams rather than index-jumps. Exclude regions with '^' prefix
targetsFile       Optional<File>           (-T) similar to -R but streams rather than index-jumps. Exclude regions with '^' prefix
threads           Optional<Integer>        number of extra output compression threads [0]
trimAltAlleles    Optional<Boolean>        (-a) trim alternate alleles not seen in the subset
noUpdate          Optional<Boolean>        (-I) do not (re)calculate INFO fields for the subset (currently INFO/AC and INFO/AN)
samples           Optional<Array<String>>  (-s) comma separated list of samples to include (or exclude with '^' prefix)
samplesFile       Optional<File>           (-S) file of samples to include (or exclude with '^' prefix)
forceSamples      Optional<Boolean>        only warn about unknown subset samples
minAc             Optional<Integer>        (-c) minimum count for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
maxAc             Optional<Integer>        (-C) maximum count for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
applyFilters      Optional<Array<String>>  (-f) require at least one of the listed FILTER strings (e.g. 'PASS,.'')
genotype          Optional<String>         (-g) [<hom|het|miss>] require one or more hom/het/missing genotype or, if prefixed with '^', exclude sites with hom/het/missing genotypes
include           Optional<String>         (-i) select sites for which the expression is true (see man page for details)
exclude           Optional<String>         (-e) exclude sites for which the expression is true (see man page for details)
known             Optional<Boolean>        (-k) select known sites only (ID is not/is '.')
novel             Optional<Boolean>        (-n) select novel sites only (ID is not/is '.')
minAlleles        Optional<Integer>        (-m) minimum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)
maxAlleles        Optional<Integer>        (-M) maximum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)
phased            Optional<Boolean>        (-p) select sites where all samples are phased
excludePhased     Optional<Boolean>        (-P) exclude sites where all samples are phased
minAf             Optional<Float>          (-q) minimum frequency for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
maxAf             Optional<Float>          (-Q) maximum frequency for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
uncalled          Optional<Boolean>        (-u) select sites without a called genotype
excludeUncalled   Optional<Boolean>        (-U) exclude sites without a called genotype
types             Optional<Array<String>>  (-v) select comma-separated list of variant types: snps,indels,mnps,other [null]
excludeTypes      Optional<Array<String>>  (-V) exclude comma-separated list of variant types: snps,indels,mnps,other [null]
private           Optional<Boolean>        (-x) select sites where the non-reference alleles are exclusive (private) to the subset samples
excludePrivate    Optional<Boolean>        (-X) exclude sites where the non-reference alleles are exclusive (private) to the subset samples
================  =======================  ==============================================================================================================================================================================

