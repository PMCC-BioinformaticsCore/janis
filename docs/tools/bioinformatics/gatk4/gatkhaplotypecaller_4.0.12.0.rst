:orphan:


GATK4: Haplotype Caller
=============================================

Description
-------------

Tool identifier: ``GatkHaplotypeCaller``

Tool path: ``janis_bioinformatics.tools.gatk4.haplotypecaller.latest import Gatk4HaplotypeCallerLatest``

Version: 4.0.12.0

Container: ``broadinstitute/gatk:4.0.12.0``



Documentation
-------------

URL
******
`https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php# <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php#>`_

Tool documentation
******************
Call germline SNPs and indels via local re-assembly of haplotypes
    
The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes 
in an active region. In other words, whenever the program encounters a region showing signs of variation, it 
discards the existing mapping information and completely reassembles the reads in that region. This allows the 
HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when 
they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at 
calling indels than position-based callers like UnifiedGenotyper.

In the GVCF workflow used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to 
generate an intermediate GVCF (not to be used in final analysis), which can then be used in GenotypeGVCFs for joint 
genotyping of multiple samples in a very efficient way. The GVCF workflow enables rapid incremental processing of 
samples as they roll off the sequencer, as well as scaling to very large cohort sizes (e.g. the 92K exomes of ExAC).

In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. 
Note however that the algorithms used to calculate variant likelihoods is not well suited to extreme allele 
frequencies (relative to ploidy) so its use is not recommended for somatic (cancer) variant discovery. 
For that purpose, use Mutect2 instead.

Finally, HaplotypeCaller is also able to correctly handle the splice junctions that make RNAseq a challenge 
for most variant callers, on the condition that the input read data has previously been processed according 
to our recommendations as documented (https://software.broadinstitute.org/gatk/documentation/article?id=4067).

Outputs
-------
======  ======  ===================================================================================================
name    type    documentation
======  ======  ===================================================================================================
out     VCFIDX  A raw, unfiltered, highly sensitive callset in VCF format. File to which variants should be written
======  ======  ===================================================================================================

Inputs
------
Find the inputs below

Required inputs
***************

=========  =============  ===========  ==========  ==================================
name       type           prefix         position  documentation
=========  =============  ===========  ==========  ==================================
inputRead  BamPair        --input                  BAM/SAM/CRAM file containing reads
reference  FastaWithDict  --reference           5  Reference sequence file
dbsnp      vcf-gz-tbi     --dbsnp               7  (Also: -D) A dbSNP VCF file.
=========  =============  ===========  ==========  ==================================

Optional inputs
***************

========================================  =======================  ===============================================  ==========  =================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                                      type                     prefix                                             position  documentation
========================================  =======================  ===============================================  ==========  =================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
activityProfileOut                        Optional<String>         --activity-profile-out                                       Output the raw activity profile results in IGV format (default: null)
alleles                                   Optional<File>           --alleles                                                    (default: null) The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES
annotateWithNumDiscoveredAlleles          Optional<Boolean>        --annotate-with-num-discovered-alleles                       If provided, we will annotate records with the number of alternate alleles that were discovered (but not necessarily genotyped) at a given site
annotation                                Optional<Array<String>>  --annotation                                                 -A: One or more specific annotations to add to variant calls
annotationGroup                           Optional<Array<String>>  --annotation-group                                           -G	One or more groups of annotations to apply to variant calls
annotationsToExclude                      Optional<Array<String>>  --annotations-to-exclude                                     -AX	One or more specific annotations to exclude from variant calls
arguments_file                            Optional<Array<File>>    --arguments_file                                             read one or more arguments files and add them to the command line
assemblyRegionOut                         Optional<String>         --assembly-region-out                                        (default: null) Output the assembly region to this IGV formatted file. Which annotations to exclude from output in the variant calls. Note that this argument has higher priority than the -A or -G arguments, so these annotations will be excluded even if they are explicitly included with the other options.
baseQualityScoreThreshold                 Optional<Integer>        --base-quality-score-threshold                               (default: 18) Base qualities below this threshold will be reduced to the minimum (6)
cloudIndexPrefetchBuffer                  Optional<Integer>        --cloud-index-prefetch-buffer                                -CIPB (default: -1) Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.
cloudPrefetchBuffer                       Optional<Integer>        --cloud-prefetch-buffer                                      -CPB (default: 40) Size of the cloud-only prefetch buffer (in MB; 0 to disable).
contaminationFractionToFilter             Optional<Double>         --contamination-fraction-to-filter                           -contamination (default: 0.0) Fraction of contamination in sequencing data (for all samples) to aggressively remove
correctOverlappingQuality                 Optional<Boolean>        --correct-overlapping-quality                                Undocumented option
disableBamIndexCaching                    Optional<Boolean>        --disable-bam-index-caching                                  -DBIC. If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified. Caching is automatically disabled if there are no intervals specified.
founderId                                 Optional<Array<String>>  --founder-id                                                 Samples representing the population "founders"
genotypingMode                            Optional<String>         --genotyping-mode                                            (default: DISCOVERY) Specifies how to determine the alternate alleles to use for genotyping. The --genotyping-mode argument is an enumerated type (GenotypingOutputMode), which can have one of the following values: DISCOVERY (The genotyper will choose the most likely alternate allele) or GENOTYPE_GIVEN_ALLELES (Only the alleles passed by the user should be considered).
heterozygosity                            Optional<Double>         --heterozygosity                                             (default: 0.001) Heterozygosity value used to compute prior likelihoods for any locus. The expected heterozygosity value used to compute prior probability that a locus is non-reference. The default priors are for provided for humans: het = 1e-3 which means that the probability of N samples being hom-ref at a site is: 1 - sum_i_2N (het / i) Note that heterozygosity as used here is the population genetics concept: http://en.wikipedia.org/wiki/Zygosity#Heterozygosity_in_population_genetics . That is, a hets value of 0.01 implies that two randomly chosen chromosomes from the population of organisms would differ from each other (one being A and the other B) at a rate of 1 in 100 bp. Note that this quantity has nothing to do with the likelihood of any given sample having a heterozygous genotype, which in the GATK is purely determined by the probability of the observed data P(D | AB) under the model that there may be a AB het genotype. The posterior probability of this AB genotype would use the het prior, but the GATK only uses this posterior probability in determining the prob. that a site is polymorphic. So changing the het parameters only increases the chance that a site will be called non-reference across all samples, but doesn't actually change the output genotype likelihoods at all, as these aren't posterior probabilities at all. The quantity that changes whether the GATK considers the possibility of a het genotype at all is the ploidy, which determines how many chromosomes each individual in the species carries.
heterozygosityStdev                       Optional<Double>         --heterozygosity-stdev                                       (default 0.01) Standard deviation of heterozygosity for SNP and indel calling.
indelHeterozygosity                       Optional<Double>         --indel-heterozygosity                                       (default: 1.25E-4) Heterozygosity for indel calling. This argument informs the prior probability of having an indel at a site. (See heterozygosity)
intervalMergingRule                       Optional<String>         --interval-merging-rule                                      -imr (default: ALL) Interval merging rule for abutting intervals. By default, the program merges abutting intervals (i.e. intervals that are directly side-by-side but do not actually overlap) into a single continuous interval. However you can change this behavior if you want them to be treated as separate intervals instead. The --interval-merging-rule argument is an enumerated type (IntervalMergingRule), which can have one of the following values:[ALL, OVERLAPPING]
maxReadsPerAlignmentStart                 Optional<Integer>        --max-reads-per-alignment-start                              (default: 50) Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.
minBaseQualityScore                       Optional<Integer>        --min-base-quality-score                                     -mbq (default: 10) Minimum base quality required to consider a base for calling
nativePairHmmThreads                      Optional<Integer>        --native-pair-hmm-threads                                    (default: 4) How many threads should a native pairHMM implementation use
nativePairHmmUseDoublePrecision           Optional<Boolean>        --native-pair-hmm-use-double-precision                       use double precision in the native pairHmm. This is slower but matches the java implementation better
numReferenceSamplesIfNoCall               Optional<Integer>        --num-reference-samples-if-no-call                           (default: 0) Number of hom-ref genotypes to infer at sites not present in a panel. When a variant is not seen in any panel, this argument controls whether to infer (and with what effective strength) that only reference alleles were observed at that site. E.g. "If not seen in 1000Genomes, treat it as AC=0, AN=2000".
outputMode                                Optional<String>         --output-mode                                                (default: EMIT_VARIANTS_ONLY) Specifies which type of calls we should output. The --output-mode argument is an enumerated type (OutputMode), which can have one of the following values: [EMIT_VARIANTS_ONLY (produces calls only at variant sites), EMIT_ALL_CONFIDENT_SITES (produces calls at variant sites and confident reference sites), EMIT_ALL_SITES (produces calls at any callable site regardless of confidence; this argument is intended only for point mutations (SNPs) in DISCOVERY mode or generally when running in GENOTYPE_GIVEN_ALLELES mode; it will by no means produce a comprehensive set of indels in DISCOVERY mode)]
pedigree                                  Optional<File>           --pedigree                                                   -ped (default: null) Pedigree file for determining the population "founders"
populationCallset                         Optional<File>           --population-callset                                         -population (default: null) Callset to use in calculating genotype priors
sampleName                                Optional<String>         --sample-name                                                -ALIAS (default: null) Name of single sample to use from a multi-sample bam. You can use this argument to specify that HC should process a single sample out of a multisample BAM file. This is especially useful if your samples are all in the same file but you need to run them individually through HC in -ERC GVC mode (which is the recommended usage). Note that the name is case-sensitive.
samplePloidy                              Optional<Integer>        --sample-ploidy                                              -ploidy (default: 2) Ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy). Sample ploidy - equivalent to number of chromosomes per pool. In pooled experiments this should be = # of samples in pool * individual sample ploidy
sitesOnlyVcfOutput                        Optional<Boolean>        --sites-only-vcf-output                                      (default: false) If true, don't emit genotype fields when writing vcf file output.
standardMinConfidenceThresholdForCalling  Optional<Double>         --standard-min-confidence-threshold-for-calling              -stand-call-conf (default: 10.0) The minimum phred-scaled confidence threshold at which variants should be called
useNewQualCalculator                      Optional<Boolean>        --use-new-qual-calculator                                    -new-qual If provided, we will use the new AF model instead of the so-called exact model
outputFilename                            Optional<Filename>       --output                                                  8  File to which variants should be written
intervals                                 Optional<bed>            --intervals                                                  -L (BASE) One or more genomic intervals over which to operate
========================================  =======================  ===============================================  ==========  =================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: Michael Franklin


*GATK4: Haplotype Caller was last updated on 2019-01-24*.
*This page was automatically generated on 2019-07-30*.
