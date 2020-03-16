:orphan:

GatkMutect2
==========================

*1 contributor · 4 versions*

USAGE: Mutect2 [arguments]
Call somatic SNVs and indels via local assembly of haplotypes
Version:4.1.2.0


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.mutect2.versions import GatkMutect2_4_1_4

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4mutect2_step",
           Gatk4Mutect2(
               tumorBams=None,
               normalBams=None,
               normalSample=None,
               reference=None,
           )
       )
       wf.output("out", source=gatk4mutect2_step.out)
   wf.output("stats", source=gatk4mutect2_step.stats)
   wf.output("f1f2r_out", source=gatk4mutect2_step.f1f2r_out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4Mutect2:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4Mutect2 > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       normalBams:
       - normalBams_0.bam
       - normalBams_1.bam
       normalSample: <value>
       reference: reference.fasta
       tumorBams:
       - tumorBams_0.bam
       - tumorBams_1.bam




5. Run Gatk4Mutect2 with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4Mutect2





Information
------------


:ID: ``Gatk4Mutect2``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.10.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.10.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.4.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24



Outputs
-----------

=========  ====================  =================
name       type                  documentation
=========  ====================  =================
out        CompressedIndexedVCF  To determine type
stats      TextFile              To determine type
f1f2r_out  CompressedTarFile     To determine type
=========  ====================  =================



Additional configuration (inputs)
---------------------------------

===================================  ==============================  ==========================================  ==========  ========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                                 type                            prefix                                        position  documentation
===================================  ==============================  ==========================================  ==========  ========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
tumorBams                            Array<IndexedBam>               -I                                                      (--input) BAM/SAM/CRAM file containing reads This argument must be specified at least once. Required.
normalBams                           Array<IndexedBam>               -I                                                      (--input) Extra BAM/SAM/CRAM file containing reads This argument must be specified at least once. Required.
normalSample                         String                          --normal-sample                                         (--normal-sample, if) May be URL-encoded as output by GetSampleName with
reference                            FastaWithIndexes                --reference                                             (-R) Reference sequence file Required.
outputFilename                       Optional<Filename>              -O                                                  20
activityProfileOut                   Optional<String>                --activity-profile-out                                  Default value: null.
addOutputSamProgramRecord            Optional<Boolean>               -add-output-sam-program-record                          (--add-output-sam-program-record)  If true, adds a PG tag to created SAM/BAM/CRAM files.  Default value: true. Possible values: {true, false}
addOutputVcfCommandLine              Optional<Boolean>               -add-output-vcf-command-line                            (--add-output-vcf-command-line)  If true, adds a command line header line to created VCF files.  Default value: true. Possible values: {true, false}
afOfAllelesNotInResource             Optional<String>                --af-of-alleles-not-in-resource                         (-default-af)  Population allele fraction assigned to alleles not found in germline resource.  Please see docs/mutect/mutect2.pdf fora derivation of the default value.  Default value: -1.0.
alleles                              Optional<String>                --alleles                                               The set of alleles for which to force genotyping regardless of evidence Default value: null.
annotation                           Optional<String>                --annotation                                            (-A) One or more specific annotations to add to variant calls This argument may be specified 0 or more times. Default value: null. Possible Values: {AlleleFraction, AS_BaseQualityRankSumTest, AS_FisherStrand, AS_InbreedingCoeff, AS_MappingQualityRankSumTest, AS_QualByDepth, AS_ReadPosRankSumTest, AS_RMSMappingQuality, AS_StrandOddsRatio, BaseQuality, BaseQualityRankSumTest, ChromosomeCounts, ClippingRankSumTest, CountNs, Coverage, DepthPerAlleleBySample, DepthPerSampleHC, ExcessHet, FisherStrand, FragmentLength, GenotypeSummaries, InbreedingCoeff, LikelihoodRankSumTest, MappingQuality, MappingQualityRankSumTest, MappingQualityZero, OrientationBiasReadCounts, OriginalAlignment, PossibleDeNovo, QualByDepth, ReadPosition, ReadPosRankSumTest, ReferenceBases, RMSMappingQuality, SampleList, StrandBiasBySample, StrandOddsRatio, TandemRepeat, UniqueAltReadCount}
annotationGroup                      Optional<String>                --annotation-group                                      (-G) One or more groups of annotations to apply to variant calls This argument may be specified 0 or more times. Default value: null. Possible Values: {AS_StandardAnnotation, ReducibleAnnotation, StandardAnnotation, StandardHCAnnotation, StandardMutectAnnotation}
annotationsToExclude                 Optional<String>                --annotations-to-exclude                                (-AX)  One or more specific annotations to exclude from variant calls  This argument may be specified 0 or more times. Default value: null. Possible Values: {BaseQuality, Coverage, DepthPerAlleleBySample, DepthPerSampleHC, FragmentLength, MappingQuality, OrientationBiasReadCounts, ReadPosition, StrandBiasBySample, TandemRepeat}
arguments_file                       Optional<File>                  --arguments_file                                        read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null.
assemblyRegionOut                    Optional<String>                --assembly-region-out                                   Output the assembly region to this IGV formatted file Default value: null.
baseQualityScoreThreshold            Optional<Integer>               --base-quality-score-threshold                          Base qualities below this threshold will be reduced to the minimum (6)  Default value: 18.
callableDepth                        Optional<Integer>               --callable-depth                                        Minimum depth to be considered callable for Mutect stats. Does not affect genotyping. Default value: 10.
cloudIndexPrefetchBuffer             Optional<Integer>               --cloud-index-prefetch-buffer                           (-CIPB)  Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.  Default value: -1.
cloudPrefetchBuffer                  Optional<Integer>               --cloud-prefetch-buffer                                 (-CPB)  Size of the cloud-only prefetch buffer (in MB; 0 to disable).  Default value: 40.
createOutputBamIndex                 Optional<Boolean>               --create-output-bam-index                               (-OBI)  If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.  Default value: true. Possible values: {true, false}
createOutputBamMd5                   Optional<Boolean>               --create-output-bam-md5                                 (-OBM)  If true, create a MD5 digest for any BAM/SAM/CRAM file created  Default value: false. Possible values: {true, false}
createOutputVariantIndex             Optional<Boolean>               --create-output-variant-index                           (-OVI)  If true, create a VCF index when writing a coordinate-sorted VCF file.  Default value: true. Possible values: {true, false}
createOutputVariantMd5               Optional<Boolean>               --create-output-variant-md5                             (-OVM)  If true, create a a MD5 digest any VCF file created.  Default value: false. Possible values: {true, false}
disableBamIndexCaching               Optional<Boolean>               --disable-bam-index-caching                             (-DBIC)  If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.  Default value: false. Possible values: {true, false}
disableReadFilter                    Optional<Boolean>               --disable-read-filter                                   (-DF)  Read filters to be disabled before analysis  This argument may be specified 0 or more times. Default value: null. Possible Values: {GoodCigarReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotSecondaryAlignmentReadFilter, PassesVendorQualityCheckReadFilter, ReadLengthReadFilter, WellformedReadFilter}
disableSequenceDictionaryValidation  Optional<Boolean>               -disable-sequence-dictionary-validation                 (--disable-sequence-dictionary-validation)  If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!  Default value: false. Possible values: {true, false}
downsamplingStride                   Optional<Integer>               --downsampling-stride                                   (-stride)  Downsample a pool of reads starting within a range of one or more bases.  Default value: 1.
excludeIntervals                     Optional<Boolean>               --exclude-intervals                                     (-XLOne) This argument may be specified 0 or more times. Default value: null.
f1r2MaxDepth                         Optional<Integer>               --f1r2-max-depth                                        sites with depth higher than this value will be grouped Default value: 200.
f1r2MedianMq                         Optional<Integer>               --f1r2-median-mq                                        skip sites with median mapping quality below this value Default value: 50.
f1r2MinBq                            Optional<Integer>               --f1r2-min-bq                                           exclude bases below this quality from pileup Default value: 20.
f1r2TarGz_outputFilename             Optional<Filename>              --f1r2-tar-gz                                           If specified, collect F1R2 counts and output files into this tar.gz file Default value: null.
founderId                            Optional<String>                -founder-id                                             (--founder-id)  Samples representing the population founders This argument may be specified 0 or more times. Default value: null.
gatkConfigFile                       Optional<String>                --gatk-config-file                                      A configuration file to use with the GATK. Default value: null.
gcsRetries                           Optional<Integer>               -gcs-retries                                            (--gcs-max-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20.
gcsProjectForRequesterPays           Optional<String>                --gcs-project-for-requester-pays                        Project to bill when accessing requester pays buckets. If unset, these buckets cannot be accessed.  Default value: .
genotypeGermlineSites                Optional<Boolean>               --genotype-germline-sites                               (EXPERIMENTAL) Call all apparent germline site even though they will ultimately be filtered.  Default value: false. Possible values: {true, false}
genotypePonSites                     Optional<Boolean>               --genotype-pon-sites                                    Call sites in the PoN even though they will ultimately be filtered. Default value: false. Possible values: {true, false}
germlineResource                     Optional<CompressedIndexedVCF>  --germline-resource                                     Population vcf of germline sequencing containing allele fractions.  Default value: null.
graph                                Optional<String>                -graph                                                  (--graph-output) Write debug assembly graph information to this file Default value: null.
help                                 Optional<Boolean>               -h                                                      (--help) display the help message Default value: false. Possible values: {true, false}
ignoreItrArtifacts                   Optional<String>                --ignore-itr-artifactsTurn                              inverted tandem repeats.  Default value: false. Possible values: {true, false}
initialTumorLod                      Optional<String>                --initial-tumor-lod                                     (-init-lod)  Log 10 odds threshold to consider pileup active.  Default value: 2.0.
intervalExclusionPadding             Optional<String>                --interval-exclusion-padding                            (-ixp)  Amount of padding (in bp) to add to each interval you are excluding.  Default value: 0.
imr                                  Optional<String>                --interval-merging-rule                                 (--interval-merging-rule)  Interval merging rule for abutting intervals  Default value: ALL. Possible values: {ALL, OVERLAPPING_ONLY}
ip                                   Optional<String>                -ipAmount                                               (--interval-padding) Default value: 0.
isr                                  Optional<String>                --interval-set-rule                                     (--interval-set-rule)  Set merging approach to use for combining interval inputs  Default value: UNION. Possible values: {UNION, INTERSECTION}
intervals                            Optional<bed>                   --intervals                                             (-L) One or more genomic intervals over which to operate This argument may be specified 0 or more times. Default value: null.
le                                   Optional<Boolean>               -LE                                                     (--lenient) Lenient processing of VCF files Default value: false. Possible values: {true, false}
maxPopulationAf                      Optional<String>                --max-population-af                                     (-max-af)  Maximum population allele frequency in tumor-only mode.  Default value: 0.01.
maxReadsPerAlignmentStart            Optional<Integer>               --max-reads-per-alignment-start                         Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.  Default value: 50.
minBaseQualityScore                  Optional<String>                --min-base-quality-score                                (-mbq:Byte)  Minimum base quality required to consider a base for calling  Default value: 10.
mitochondriaMode                     Optional<Boolean>               --mitochondria-mode                                     Mitochondria mode sets emission and initial LODs to 0. Default value: false. Possible values: {true, false}
nativePairHmmThreads                 Optional<Integer>               --native-pair-hmm-threads                               How many threads should a native pairHMM implementation use  Default value: 4.
nativePairHmmUseDoublePrecision      Optional<Boolean>               --native-pair-hmm-use-double-precision                  use double precision in the native pairHmm. This is slower but matches the java implementation better  Default value: false. Possible values: {true, false}
normalLod                            Optional<Double>                --normal-lod                                            Log 10 odds threshold for calling normal variant non-germline. Default value: 2.2.
encode                               Optional<String>                -encode                                                 This argument may be specified 0 or more times. Default value: null.
panelOfNormals                       Optional<CompressedIndexedVCF>  --panel-of-normals                                      (--panel-of-normals)  VCF file of sites observed in normal.  Default value: null.
pcrIndelQual                         Optional<Integer>               --pcr-indel-qual                                        Phred-scaled PCR SNV qual for overlapping fragments Default value: 40.
pcrSnvQual                           Optional<Integer>               --pcr-snv-qual                                          Phred-scaled PCR SNV qual for overlapping fragments Default value: 40.
pedigree                             Optional<String>                --pedigree                                              (-ped) Pedigree file for determining the population founders. Default value: null.
quiet                                Optional<Boolean>               --QUIET                                                 Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false}
readFilter                           Optional<String>                --read-filter                                           (-RF) Read filters to be applied before analysis This argument may be specified 0 or more times. Default value: null. Possible Values: {AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter}
readIndex                            Optional<String>                -read-index                                             (--read-index)  Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically.  This argument may be specified 0 or more times. Default value: null.
readValidationStringency             Optional<String>                --read-validation-stringency                            (-VS:ValidationStringency)  Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: SILENT. Possible values: {STRICT, LENIENT, SILENT}
secondsBetweenProgressUpdates        Optional<Double>                -seconds-between-progress-updates                       (--seconds-between-progress-updates)  Output traversal statistics every time this many seconds elapse  Default value: 10.0.
sequenceDictionary                   Optional<String>                -sequence-dictionary                                    (--sequence-dictionary)  Use the given sequence dictionary as the master/canonical sequence dictionary.  Must be a .dict file.  Default value: null.
sitesOnlyVcfOutput                   Optional<Boolean>               --sites-only-vcf-output                                 If true, don't emit genotype fields when writing vcf file output.  Default value: false. Possible values: {true, false}
tmpDir                               Optional<String>                --tmp-dir                                               Temp directory to use. Default value: null.
tumorLodToEmit                       Optional<String>                --tumor-lod-to-emit                                     (-emit-lod)  Log 10 odds threshold to emit variant to VCF.  Default value: 3.0.
tumor                                Optional<String>                -tumor                                                  (--tumor-sample) BAM sample name of tumor. May be URL-encoded as output by GetSampleName with -encode argument.  Default value: null.
jdkDeflater                          Optional<Boolean>               -jdk-deflater                                           (--use-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false}
jdkInflater                          Optional<Boolean>               -jdk-inflater                                           (--use-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false}
verbosity                            Optional<String>                -verbosity                                              (--verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG}
version                              Optional<Boolean>               --version                                               display the version number for this tool Default value: false. Possible values: {true, false}
activeProbabilityThreshold           Optional<Double>                --active-probability-threshold                          Minimum probability for a locus to be considered active.  Default value: 0.002.
adaptivePruningInitialErrorRate      Optional<Double>                --adaptive-pruning-initial-error-rate                   Initial base error rate estimate for adaptive pruning  Default value: 0.001.
allowNonUniqueKmersInRef             Optional<Boolean>               --allow-non-unique-kmers-in-ref                         Allow graphs that have non-unique kmers in the reference  Default value: false. Possible values: {true, false}
assemblyRegionPadding                Optional<Integer>               --assembly-region-padding                               Number of additional bases of context to include around each assembly region  Default value: 100.
bamout                               Optional<String>                -bamout                                                 (--bam-output) File to which assembled haplotypes should be written Default value: null.
bamWriterType                        Optional<String>                --bam-writer-type                                       Which haplotypes should be written to the BAM Default value: CALLED_HAPLOTYPES. Possible values: {ALL_POSSIBLE_HAPLOTYPES, CALLED_HAPLOTYPES}
debugAssembly                        Optional<String>                --debug-assembly                                        (-debug)  Print out verbose debug information about each assembly region  Default value: false. Possible values: {true, false}
disableAdaptivePruning               Optional<Boolean>               --disable-adaptive-pruning                              Disable the adaptive algorithm for pruning paths in the graph  Default value: false. Possible values: {true, false}
disableToolDefaultAnnotations        Optional<Boolean>               -disable-tool-default-annotations                       (--disable-tool-default-annotations)  Disable all tool default annotations  Default value: false. Possible values: {true, false}
disableToolDefaultReadFilters        Optional<Boolean>               -disable-tool-default-read-filters                      (--disable-tool-default-read-filters)  Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)  Default value: false. Possible values: {true, false}
dontIncreaseKmerSizesForCycles       Optional<Boolean>               --dont-increase-kmer-sizes-for-cycles                   Disable iterating over kmer sizes when graph cycles are detected  Default value: false. Possible values: {true, false}
dontTrimActiveRegions                Optional<Boolean>               --dont-trim-active-regions                              If specified, we will not trim down the active region from the full region (active + extension) to just the active interval for genotyping  Default value: false. Possible values: {true, false}
dontUseSoftClippedBases              Optional<Boolean>               --dont-use-soft-clipped-bases                           Do not analyze soft clipped bases in the reads  Default value: false. Possible values: {true, false}
erc                                  Optional<String>                -ERC                                                    (--emit-ref-confidence)  (BETA feature) Mode for emitting reference confidence scores  Default value: NONE. Possible values: {NONE, BP_RESOLUTION, GVCF}
enableAllAnnotations                 Optional<Boolean>               --enable-all-annotations                                Use all possible annotations (not for the faint of heart)  Default value: false. Possible values: {true, false}
forceActive                          Optional<Boolean>               --force-active                                          If provided, all regions will be marked as active Default value: false. Possible values: {true, false}
genotypeFilteredAlleles              Optional<Boolean>               --genotype-filtered-alleles                             Whether to force genotype even filtered alleles  Default value: false. Possible values: {true, false}
gvcfLodBand                          Optional<String>                --gvcf-lod-band                                         (-LODB) Exclusive upper bounds for reference confidence LOD bands (must be specified in increasing order)  This argument may be specified 0 or more times. Default value: [-2.5, -2.0, -1.5,
kmerSize                             Optional<Integer>               --kmer-size                                             Kmer size to use in the read threading assembler This argument may be specified 0 or more times. Default value: [10, 25].
maxAssemblyRegionSize                Optional<Integer>               --max-assembly-region-size                              Maximum size of an assembly region  Default value: 300.
mnpDist                              Optional<Integer>               -mnp-dist                                               (--max-mnp-distance)  Two or more phased substitutions separated by this distance or less are merged into MNPs.  Default value: 1.
maxNumHaplotypesInPopulation         Optional<Integer>               --max-num-haplotypes-in-population                      Maximum number of haplotypes to consider for your population  Default value: 128.
maxProbPropagationDistance           Optional<Integer>               --max-prob-propagation-distance                         Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions  Default value: 50.
maxSuspiciousReadsPerAlignmentStart  Optional<Integer>               --max-suspicious-reads-per-alignment-start              Maximum number of suspicious reads (mediocre mapping quality or too many substitutions) allowed in a downsampling stride.  Set to 0 to disable.  Default value: 0.
maxUnprunedVariants                  Optional<Integer>               --max-unpruned-variants                                 Maximum number of variants in graph the adaptive pruner will allow  Default value: 100.
minAssemblyRegionSize                Optional<Integer>               --min-assembly-region-size                              Minimum size of an assembly region  Default value: 50.
minDanglingBranchLength              Optional<Integer>               --min-dangling-branch-length                            Minimum length of a dangling branch to attempt recovery  Default value: 4.
minPruning                           Optional<Integer>               --min-pruning                                           Minimum support to not prune paths in the graph Default value: 2.
minimumAlleleFraction                Optional<Float>                 --minimum-allele-fraction                               (-min-AF)  Lower bound of variant allele fractions to consider when calculating variant LOD  Default value: 0.0.
numPruningSamples                    Optional<Integer>               --num-pruning-samples                                   Default value: 1.
pairHmmGapContinuationPenalty        Optional<Integer>               --pair-hmm-gap-continuation-penalty                     Flat gap continuation penalty for use in the Pair HMM  Default value: 10.
pairhmm                              Optional<String>                -pairHMM                                                (--pair-hmm-implementation)  The PairHMM implementation to use for genotype likelihood calculations  Default value: FASTEST_AVAILABLE. Possible values: {EXACT, ORIGINAL, LOGLESS_CACHING, AVX_LOGLESS_CACHING, AVX_LOGLESS_CACHING_OMP, EXPERIMENTAL_FPGA_LOGLESS_CACHING, FASTEST_AVAILABLE}
pcrIndelModel                        Optional<String>                --pcr-indel-model                                       The PCR indel model to use  Default value: CONSERVATIVE. Possible values: {NONE, HOSTILE, AGGRESSIVE, CONSERVATIVE}
phredScaledGlobalReadMismappingRate  Optional<Integer>               --phred-scaled-global-read-mismapping-rate              The global assumed mismapping rate for reads  Default value: 45.
pruningLodThreshold                  Optional<Float>                 --pruning-lod-thresholdLn                               Default value: 2.302585092994046.
recoverAllDanglingBranches           Optional<Boolean>               --recover-all-dangling-branches                         Recover all dangling branches  Default value: false. Possible values: {true, false}
showhidden                           Optional<Boolean>               -showHidden                                             (--showHidden)  display hidden arguments  Default value: false. Possible values: {true, false}
smithWaterman                        Optional<String>                --smith-waterman                                        Which Smith-Waterman implementation to use, generally FASTEST_AVAILABLE is the right choice  Default value: JAVA. Possible values: {FASTEST_AVAILABLE, AVX_ENABLED, JAVA}
ambigFilterBases                     Optional<Integer>               --ambig-filter-bases                                    Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction.  Default value: null.  Cannot be used in conjuction with argument(s) maxAmbiguousBaseFraction
ambigFilterFrac                      Optional<Double>                --ambig-filter-frac                                     Threshold fraction of ambiguous bases Default value: 0.05. Cannot be used in conjuction with argument(s) maxAmbiguousBases
maxFragmentLength                    Optional<Integer>               --max-fragment-length                                   Default value: 1000000.
minFragmentLength                    Optional<Integer>               --min-fragment-length                                   Default value: 0.
keepIntervals                        Optional<String>                --keep-intervals                                        One or more genomic intervals to keep This argument must be specified at least once. Required.
library                              Optional<String>                -library                                                (--library) Name of the library to keep This argument must be specified at least once. Required.
maximumMappingQuality                Optional<Integer>               --maximum-mapping-quality                               Maximum mapping quality to keep (inclusive)  Default value: null.
minimumMappingQuality                Optional<Integer>               --minimum-mapping-quality                               Minimum mapping quality to keep (inclusive)  Default value: 20.
dontRequireSoftClipsBothEnds         Optional<Boolean>               --dont-require-soft-clips-both-ends                     Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block  Default value: false. Possible values: {true, false}
filterTooShort                       Optional<Integer>               --filter-too-short                                      Minimum number of aligned bases Default value: 30.
platformFilterName                   Optional<String>                --platform-filter-name                                  This argument must be specified at least once. Required.
blackListedLanes                     Optional<String>                --black-listed-lanes                                    Platform unit (PU) to filter out This argument must be specified at least once. Required.
readGroupBlackList                   Optional<String>                --read-group-black-listThe                              This argument must be specified at least once. Required.
keepReadGroup                        Optional<String>                --keep-read-group                                       The name of the read group to keep Required.
maxReadLength                        Optional<Integer>               --max-read-length                                       Keep only reads with length at most equal to the specified value Default value: 2147483647.
minReadLength                        Optional<Integer>               --min-read-length                                       Keep only reads with length at least equal to the specified value Default value: 30.
readName                             Optional<String>                --read-name                                             Keep only reads with this read name Required.
keepReverseStrandOnly                Optional<Boolean>               --keep-reverse-strand-only                              Keep only reads on the reverse strand  Required. Possible values: {true, false}
sample                               Optional<String>                -sample                                                 (--sample) The name of the sample(s) to keep, filtering out all others This argument must be specified at least once. Required.
===================================  ==============================  ==========================================  ==========  ========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
