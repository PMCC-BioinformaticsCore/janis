:orphan:

freebayes
=========

*2 contributors · 2 versions*

usage: freebayes [OPTION] ... [BAM FILE] ...
Bayesian haplotype-based polymorphism discovery.
Version:1.3.1


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.freebayes.versions import FreeBayes_1_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "freebayes_step",
           freebayes(
               bams=None,
               reference=None,
               theta=None,
               ploidy=None,
               refQual=None,
               maxNumOfAlleles=None,
               haplotypeLength=None,
               minRepSize=None,
               minRepEntropy=None,
               useDupFlag=None,
               minMappingQual=None,
               minBaseQual=None,
               minSupQsum=None,
               minSupMQsum=None,
               minSupBQthres=None,
               maxMisMatchFrac=None,
               minAltFrac=None,
               minAltCount=None,
               minAltQSum=None,
               minAltTotal=None,
               minCov=None,
               probContamin=None,
               genotypingMaxIter=None,
               genotypingMaxBDepth=None,
               postIntegrationLim=None,
               readDepFact=None,
           )
       )
       wf.output("out", source=freebayes_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for freebayes:

.. code-block:: bash

   # user inputs
   janis inputs freebayes > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bams:
       - bams_0.bam
       - bams_1.bam
       maxNumOfAlleles: 0
       minAltQSum: 0
       minBaseQual: 0
       minCov: 0
       minSupMQsum: 0
       minSupQsum: 0
       reference: reference.fasta
       useDupFlag: true




5. Run freebayes with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       freebayes





Information
------------


:ID: ``freebayes``
:URL: `https://github.com/ekg/freebayes <https://github.com/ekg/freebayes>`_
:Versions: 1.3.1, 1.2
:Container: shollizeck/freebayes:1.3.1
:Authors: Sebastian Hollizeck, Michael Franklin
:Citations: Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907 [q-bio.GN] 2012
:Created: 2019-10-19
:Updated: 2019-10-19



Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============



Additional configuration (inputs)
---------------------------------

=====================  ==================  ================================  ==========  =============================================================================================================================================================================================================================================================================================================
name                   type                prefix                            position    documentation
=====================  ==================  ================================  ==========  =============================================================================================================================================================================================================================================================================================================
bams                   Array<IndexedBam>   -b                                            Add FILE to the set of BAM files to be analyzed.
reference              FastaWithIndexes    -f                                            Use FILE as the reference sequence for analysis. An index file (FILE.fai) will be created if none exists. If neither --targets nor --region are specified, FreeBayes will analyze every position in this reference.
theta                  Float               -T                                            The expected mutation rate or pairwise nucleotide diversity among the population under analysis. This serves as the single parameter to the Ewens Sampling Formula prior model default: 0.001
ploidy                 Integer             -p                                            Sets the default ploidy for the analysis to N. default: 2
refQual                String              --reference-quality                           --reference-quality MQ,BQ  Assign mapping quality of MQ to the reference allele at each site and base quality of BQ. default: 100,60
maxNumOfAlleles        Integer             -n                                            Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores. (Set to 0 to use all; default: all)
haplotypeLength        Integer             --haplotype-length                            Allow haplotype calls with contiguous embedded matches of up to this length. Set N=-1 to disable clumping. (default: 3)
minRepSize             Integer             --min-repeat-size                             When assembling observations across repeats, require the total repeat length at least this many bp. (default: 5)
minRepEntropy          Integer             --min-repeat-entropy                          To detect interrupted repeats, build across sequence until it has  entropy > N bits per bp. Set to 0 to turn off. (default: 1)
useDupFlag             Boolean             -4                                            Include duplicate-marked alignments in the analysis. default: exclude duplicates marked as such in alignments
minMappingQual         Integer             -m                                            Exclude alignments from analysis if they have a mapping quality less than Q. default: 1
minBaseQual            Integer             -q                                            -q --min-base-quality Q Exclude alleles from analysis if their supporting base quality is less than Q. default: 0
minSupQsum             Integer             -R                                            -R --min-supporting-allele-qsum Q Consider any allele in which the sum of qualities of supporting observations is at least Q. default: 0
minSupMQsum            Integer             -Y                                            -Y --min-supporting-mapping-qsum Q Consider any allele in which and the sum of mapping qualities of supporting reads is at least Q. default: 0
minSupBQthres          Integer             -Q                                            -Q --mismatch-base-quality-threshold Q Count mismatches toward --read-mismatch-limit if the base quality of the mismatch is >= Q. default: 10
maxMisMatchFrac        Float               -z                                            -z --read-max-mismatch-fraction N Exclude reads with more than N [0,1] fraction of mismatches where each mismatch has base quality >= mismatch-base-quality-threshold default: 1.0
minAltFrac             Float               -F                                            -F --min-alternate-fraction N Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position. default: 0.05
minAltCount            Integer             -C                                            -C --min-alternate-count N Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position. default: 2
minAltQSum             Integer             -3                                            -3 --min-alternate-qsum N Require at least this sum of quality of observations supporting an alternate allele within a single individual in order to evaluate the position. default: 0
minAltTotal            Integer             -G                                            -G --min-alternate-total N Require at least this count of observations supporting an alternate allele within the total population in order to use the allele in analysis. default: 1
minCov                 Integer             --min-coverage                                --min-coverage N Require at least this coverage to process a site. default: 0
probContamin           Float               --prob-contamination                          --prob-contamination F An estimate of contamination to use for all samples. default: 10e-9
genotypingMaxIter      Integer             -B                                            -B --genotyping-max-iterations N Iterate no more than N times during genotyping step. default: 1000.
genotypingMaxBDepth    Integer             --genotyping-max-banddepth                    --genotyping-max-banddepth N Integrate no deeper than the Nth best genotype by likelihood when genotyping. default: 6.
postIntegrationLim     String              -W                                            -W --posterior-integration-limits N,M Integrate all genotype combinations in our posterior space which include no more than N samples with their Mth best data likelihood. default: 1,3.
readDepFact            Float               -D                                            -D --read-dependence-factor N Incorporate non-independence of reads by scaling successive observations by this factor during data likelihood calculations. default: 0.9
bamList                Optional<TextFile>  -L                                            A file containing a list of BAM files to be analyzed.
targetsFile            Optional<bed>       -t                                            Limit analysis to targets listed in the BED-format FILE.
region                 Optional<String>    -r                                            <chrom>:<start_position>-<end_position> Limit analysis to the specified region, 0-base coordinates, end_position not included (same as BED format). Either '-' or '..' maybe used as a separator.
samplesFile            Optional<TextFile>  -s                                            FILE  Limit analysis to samples listed (one per line) in the FILE. By default FreeBayes will analyze all samples in its input BAM files.
popFile                Optional<TextFile>  --populations                                 FILE Each line of FILE should list a sample and a population which it is part of. The population-based bayesian inference model will then be partitioned on the basis of the populations.
cnvFile                Optional<TextFile>  -A                                            FILE Read a copy number map from the BED file FILE, which has either a sample-level ploidy: sample name, copy number or a region-specific format: reference sequence, start, end, sample name, copy number ... for each region in each sample which does not have the default copy number as set by --ploidy.
outputFilename         Optional<Filename>  -v                                            FILE Output VCF-format results to FILE. (default: stdout)
gvcfFlag               Optional<Boolean>   --gvcf                                        Write gVCF output, which indicates coverage in uncalled regions.
gvcfChunkSize          Optional<Integer>   --gvcf-chunk                                  When writing gVCF output emit a record for every NUM bases.
candidateVcf           Optional<File>      -@                                            Use variants reported in VCF file as input to the algorithm. Variants in this file will included in the output even if there is not enough support in the data to pass input filters.
restrictSitesFlag      Optional<Boolean>   -l                                            Only provide variant calls and genotype likelihoods for sites and alleles which are provided in the VCF input, and provide output in the VCF for all input alleles, not just those which have support in the data.
candidateHaploVcf      Optional<File>      --haplotype-basis-alleles                     When specified, only variant alleles provided in this input VCF will be used for the construction of complex or haplotype alleles.
reportHapAllelesFlag   Optional<Boolean>   --report-all-haplotype-alleles                At sites where genotypes are made over haplotype alleles, provide information about all alleles in output, not only those which are called.
monomorphicFlag        Optional<Boolean>   --report-monomorphic                          Report even loci which appear to be monomorphic, and report all considered alleles, even those which are not in called genotypes. Loci which do not have any potential alternates have '.' for ALT.
polyMoprhProbFlag      Optional<Float>     -P                                            Report sites if the probability that there is a polymorphism at the site is greater than N. default: 0.0. Note that post-filtering is generally recommended over the use of this parameter.
strictFlag             Optional<Boolean>   --strict-vcf                                  Generate strict VCF format (FORMAT/GQ will be an int)
pooledDiscreteFlag     Optional<Boolean>   -J                                            Assume that samples result from pooled sequencing. Model pooled samples using discrete genotypes across pools. When using this flag, set --ploidy to the number of alleles in each sample or use the --cnv-map to define per-sample ploidy.
pooledContinousFlag    Optional<Boolean>   -K                                            Output all alleles which pass input filters, regardles of genotyping outcome or model.
addRefFlag             Optional<Boolean>   -Z                                            This flag includes the reference allele in the analysis as if it is another sample from the same population.
ignoreSNPsFlag         Optional<Boolean>   -I                                            Ignore SNP alleles.
ignoreINDELsFlag       Optional<Boolean>   -i                                            Ignore insertion and deletion alleles.
ignoreMNPsFlag         Optional<Boolean>   -X                                            Ignore multi-nuceotide polymorphisms, MNPs.
ignoreComplexVarsFlag  Optional<Boolean>   -u                                            Ignore complex events (composites of other classes).
maxNumOfComplexVars    Optional<Integer>   -E
noPartObsFlag          Optional<Boolean>   --no-partial-observations                     Exclude observations which do not fully span the dynamically-determined detection window. (default, use all observations, dividing partial support across matching haplotypes when generating haplotypes.)
noNormaliseFlag        Optional<Boolean>   -O                                            Turn off left-alignment of indels, which is enabled by default.
readMisMatchLim        Optional<Integer>   -U                                            -U --read-mismatch-limit N Exclude reads with more than N mismatches where each mismatch has base quality >= mismatch-base-quality-threshold. default: ~unbounded
readSNPLim             Optional<Integer>   -$                                            -$ --read-snp-limit N Exclude reads with more than N base mismatches, ignoring gaps with quality >= mismatch-base-quality-threshold. default: ~unbounded
readINDELLim           Optional<Integer>   -e                                            -e --read-indel-limit N Exclude reads with more than N separate gaps. default: ~unbounded
standardFilterFlag     Optional<Boolean>   -0                                            -0 --standard-filters Use stringent input base and mapping quality filters Equivalent to -m 30 -q 20 -R 0 -S 0
maxCov                 Optional<Integer>   --limit-coverage                              Downsample per-sample coverage to this level if greater than this coverage. default: no limit
noPopPriorsFlag        Optional<Boolean>   -k                                            -k --no-population-priors Equivalent to --pooled-discrete --hwe-priors-off and removal of Ewens Sampling Formula component of priors.
noHWEPriorsFlag        Optional<Boolean>   -w                                            -w --hwe-priors-off Disable estimation of the probability of the combination arising under HWE given the allele frequency as estimated by observation frequency.
noBinOBSPriorsFlag     Optional<Boolean>   -V                                            -V --binomial-obs-priors-off Disable incorporation of prior expectations about observations. Uses read placement probability, strand balance probability, and read position (5'-3') probability.
noABPriorsFlag         Optional<Boolean>   -a                                            -a --allele-balance-priors-off Disable use of aggregate probability of observation balance between alleles as a component of the priors.
obsBiasFile            Optional<TextFile>  --observation-bias                            --observation-bias FILE Read length-dependent allele observation biases from FILE. The format is [length] [alignment efficiency relative to reference] where the efficiency is 1 if there is no relative observation bias.
baseQualCap            Optional<Integer>   --base-quality-cap                            --base-quality-cap Q Limit estimated observation quality by capping base quality at Q.
legGLScalc             Optional<Boolean>   --legacy-gls                                  --legacy-gls Use legacy (polybayes equivalent) genotype likelihood calculations
contaminEst            Optional<TextFile>  --contamination-estimates                     --contamination-estimates FILE A file containing per-sample estimates of contamination, such as those generated by VerifyBamID. The format should be: sample p(read=R|genotype=AR) p(read=A|genotype=AA) Sample '*' can be used to set default contamination estimates.
repoprtMaxGLFlag       Optional<Boolean>   --report-genotype-likelihood-max              --report-genotype-likelihood-max Report genotypes using the maximum-likelihood estimate provided from genotype likelihoods.
excludeUnObsGT         Optional<Boolean>   -N                                            -N --exclude-unobserved-genotypes Skip sample genotypings for which the sample has no supporting reads.
gtVarThres             Optional<Integer>   -S                                            -S --genotype-variant-threshold N Limit posterior integration to samples where the second-best genotype likelihood is no more than log(N) from the highest genotype likelihood for the sample. default: ~unbounded
useMQFlag              Optional<Boolean>   -j                                            -j --use-mapping-quality Use mapping quality of alleles when calculating data likelihoods.
harmIndelQualFlag      Optional<Boolean>   -H                                            -H --harmonic-indel-quality Use a weighted sum of base qualities around an indel, scaled by the distance from the indel. By default use a minimum BQ in flanking sequence.
gtQuals                Optional<Boolean>   -=                                            -= --genotype-qualities Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output.
skipCov                Optional<Integer>   --skip-coverage                               Skip processing of alignments overlapping positions with coverage >N. This filters sites above this coverage, but will also reduce data nearby. default: no limit
=====================  ==================  ================================  ==========  =============================================================================================================================================================================================================================================================================================================
