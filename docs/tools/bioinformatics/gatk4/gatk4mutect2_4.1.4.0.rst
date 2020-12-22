:orphan:

GatkMutect2
==========================

``Gatk4Mutect2`` · *1 contributor · 7 versions*

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
           GatkMutect2_4_1_4(
               tumorBams=None,
               reference=None,
           )
       )
       wf.output("out", source=gatk4mutect2_step.out)
       wf.output("stats", source=gatk4mutect2_step.stats)
       wf.output("f1f2r_out", source=gatk4mutect2_step.f1f2r_out)
       wf.output("bam", source=gatk4mutect2_step.bam)
    

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
:Versions: 4.1.8.1, 4.1.7.0, 4.1.6.0, 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.4.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24


Outputs
-----------

=========  ====================  ====================================================
name       type                  documentation
=========  ====================  ====================================================
out        Gzipped<VCF>          To determine type
stats      TextFile              To determine type
f1f2r_out  CompressedTarFile     To determine type
bam        Optional<IndexedBam>  File to which assembled haplotypes should be written
=========  ====================  ====================================================


Additional configuration (inputs)
---------------------------------

===================================  ===========================  ==========================================  ==========  ========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                                 type                         prefix                                        position  documentation
===================================  ===========================  ==========================================  ==========  ========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
tumorBams                            Array<IndexedBam>            -I                                                      (--input) BAM/SAM/CRAM file containing reads This argument must be specified at least once. Required.
reference                            FastaWithIndexes             --reference                                             (-R) Reference sequence file Required.
javaOptions                          Optional<Array<String>>
compression_level                    Optional<Integer>                                                                    Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
normalBams                           Optional<Array<IndexedBam>>  -I                                                      (--input) Extra BAM/SAM/CRAM file containing reads This argument must be specified at least once. Required.
normalSample                         Optional<String>             --normal-sample                                         (--normal-sample, if) May be URL-encoded as output by GetSampleName with
outputPrefix                         Optional<String>                                                                     Used as a prefix for the outputFilename if not specified, with format: {outputPrefix}.vcf.gz
outputFilename                       Optional<Filename>           -O                                                  20
outputBamName                        Optional<String>             -bamout                                                 File to which assembled haplotypes should be written
activityProfileOut                   Optional<String>             --activity-profile-out                                  Default value: null.
addOutputSamProgramRecord            Optional<Boolean>            -add-output-sam-program-record                          (--add-output-sam-program-record)  If true, adds a PG tag to created SAM/BAM/CRAM files.  Default value: true. Possible values: {true, false}
addOutputVcfCommandLine              Optional<Boolean>            -add-output-vcf-command-line                            (--add-output-vcf-command-line)  If true, adds a command line header line to created VCF files.  Default value: true. Possible values: {true, false}
afOfAllelesNotInResource             Optional<String>             --af-of-alleles-not-in-resource                         (-default-af)  Population allele fraction assigned to alleles not found in germline resource.  Please see docs/mutect/mutect2.pdf fora derivation of the default value.  Default value: -1.0.
alleles                              Optional<String>             --alleles                                               The set of alleles for which to force genotyping regardless of evidence Default value: null.
annotation                           Optional<String>             --annotation                                            (-A) One or more specific annotations to add to variant calls This argument may be specified 0 or more times. Default value: null. Possible Values: {AlleleFraction, AS_BaseQualityRankSumTest, AS_FisherStrand, AS_InbreedingCoeff, AS_MappingQualityRankSumTest, AS_QualByDepth, AS_ReadPosRankSumTest, AS_RMSMappingQuality, AS_StrandOddsRatio, BaseQuality, BaseQualityRankSumTest, ChromosomeCounts, ClippingRankSumTest, CountNs, Coverage, DepthPerAlleleBySample, DepthPerSampleHC, ExcessHet, FisherStrand, FragmentLength, GenotypeSummaries, InbreedingCoeff, LikelihoodRankSumTest, MappingQuality, MappingQualityRankSumTest, MappingQualityZero, OrientationBiasReadCounts, OriginalAlignment, PossibleDeNovo, QualByDepth, ReadPosition, ReadPosRankSumTest, ReferenceBases, RMSMappingQuality, SampleList, StrandBiasBySample, StrandOddsRatio, TandemRepeat, UniqueAltReadCount}
annotationGroup                      Optional<String>             --annotation-group                                      (-G) One or more groups of annotations to apply to variant calls This argument may be specified 0 or more times. Default value: null. Possible Values: {AS_StandardAnnotation, ReducibleAnnotation, StandardAnnotation, StandardHCAnnotation, StandardMutectAnnotation}
annotationsToExclude                 Optional<String>             --annotations-to-exclude                                (-AX)  One or more specific annotations to exclude from variant calls  This argument may be specified 0 or more times. Default value: null. Possible Values: {BaseQuality, Coverage, DepthPerAlleleBySample, DepthPerSampleHC, FragmentLength, MappingQuality, OrientationBiasReadCounts, ReadPosition, StrandBiasBySample, TandemRepeat}
arguments_file                       Optional<File>               --arguments_file                                        read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null.
assemblyRegionOut                    Optional<String>             --assembly-region-out                                   Output the assembly region to this IGV formatted file Default value: null.
baseQualityScoreThreshold            Optional<Integer>            --base-quality-score-threshold                          Base qualities below this threshold will be reduced to the minimum (6)  Default value: 18.
callableDepth                        Optional<Integer>            --callable-depth                                        Minimum depth to be considered callable for Mutect stats. Does not affect genotyping. Default value: 10.
cloudIndexPrefetchBuffer             Optional<Integer>            --cloud-index-prefetch-buffer                           (-CIPB)  Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.  Default value: -1.
cloudPrefetchBuffer                  Optional<Integer>            --cloud-prefetch-buffer                                 (-CPB)  Size of the cloud-only prefetch buffer (in MB; 0 to disable).  Default value: 40.
createOutputBamIndex                 Optional<Boolean>            --create-output-bam-index                               (-OBI)  If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.  Default value: true. Possible values: {true, false}
createOutputBamMd5                   Optional<Boolean>            --create-output-bam-md5                                 (-OBM)  If true, create a MD5 digest for any BAM/SAM/CRAM file created  Default value: false. Possible values: {true, false}
createOutputVariantIndex             Optional<Boolean>            --create-output-variant-index                           (-OVI)  If true, create a VCF index when writing a coordinate-sorted VCF file.  Default value: true. Possible values: {true, false}
createOutputVariantMd5               Optional<Boolean>            --create-output-variant-md5                             (-OVM)  If true, create a a MD5 digest any VCF file created.  Default value: false. Possible values: {true, false}
disableBamIndexCaching               Optional<Boolean>            --disable-bam-index-caching                             (-DBIC)  If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.  Default value: false. Possible values: {true, false}
disableReadFilter                    Optional<Boolean>            --disable-read-filter                                   (-DF)  Read filters to be disabled before analysis  This argument may be specified 0 or more times. Default value: null. Possible Values: {GoodCigarReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotSecondaryAlignmentReadFilter, PassesVendorQualityCheckReadFilter, ReadLengthReadFilter, WellformedReadFilter}
disableSequenceDictionaryValidation  Optional<Boolean>            -disable-sequence-dictionary-validation                 (--disable-sequence-dictionary-validation)  If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!  Default value: false. Possible values: {true, false}
downsamplingStride                   Optional<Integer>            --downsampling-stride                                   (-stride)  Downsample a pool of reads starting within a range of one or more bases.  Default value: 1.
excludeIntervals                     Optional<Boolean>            --exclude-intervals                                     (-XLOne) This argument may be specified 0 or more times. Default value: null.
f1r2MaxDepth                         Optional<Integer>            --f1r2-max-depth                                        sites with depth higher than this value will be grouped Default value: 200.
f1r2MedianMq                         Optional<Integer>            --f1r2-median-mq                                        skip sites with median mapping quality below this value Default value: 50.
f1r2MinBq                            Optional<Integer>            --f1r2-min-bq                                           exclude bases below this quality from pileup Default value: 20.
f1r2TarGz_outputFilename             Optional<Filename>           --f1r2-tar-gz                                           If specified, collect F1R2 counts and output files into this tar.gz file Default value: null.
founderId                            Optional<String>             -founder-id                                             (--founder-id)  Samples representing the population founders This argument may be specified 0 or more times. Default value: null.
gatkConfigFile                       Optional<String>             --gatk-config-file                                      A configuration file to use with the GATK. Default value: null.
gcsRetries                           Optional<Integer>            -gcs-retries                                            (--gcs-max-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20.
gcsProjectForRequesterPays           Optional<String>             --gcs-project-for-requester-pays                        Project to bill when accessing requester pays buckets. If unset, these buckets cannot be accessed.  Default value: .
genotypeGermlineSites                Optional<Boolean>            --genotype-germline-sites                               (EXPERIMENTAL) Call all apparent germline site even though they will ultimately be filtered.  Default value: false. Possible values: {true, false}
genotypePonSites                     Optional<Boolean>            --genotype-pon-sites                                    Call sites in the PoN even though they will ultimately be filtered. Default value: false. Possible values: {true, false}
germlineResource                     Optional<Gzipped<VCF>>       --germline-resource                                     Population vcf of germline sequencing containing allele fractions.  Default value: null.
graph                                Optional<String>             -graph                                                  (--graph-output) Write debug assembly graph information to this file Default value: null.
help                                 Optional<Boolean>            -h                                                      (--help) display the help message Default value: false. Possible values: {true, false}
ignoreItrArtifacts                   Optional<String>             --ignore-itr-artifactsTurn                              inverted tandem repeats.  Default value: false. Possible values: {true, false}
initialTumorLod                      Optional<String>             --initial-tumor-lod                                     (-init-lod)  Log 10 odds threshold to consider pileup active.  Default value: 2.0.
intervalExclusionPadding             Optional<String>             --interval-exclusion-padding                            (-ixp)  Amount of padding (in bp) to add to each interval you are excluding.  Default value: 0.
imr                                  Optional<String>             --interval-merging-rule                                 (--interval-merging-rule)  Interval merging rule for abutting intervals  Default value: ALL. Possible values: {ALL, OVERLAPPING_ONLY}
ip                                   Optional<String>             -ipAmount                                               (--interval-padding) Default value: 0.
isr                                  Optional<String>             --interval-set-rule                                     (--interval-set-rule)  Set merging approach to use for combining interval inputs  Default value: UNION. Possible values: {UNION, INTERSECTION}
intervals                            Optional<bed>                --intervals                                             (-L) One or more genomic intervals over which to operate This argument may be specified 0 or more times. Default value: null.
le                                   Optional<Boolean>            -LE                                                     (--lenient) Lenient processing of VCF files Default value: false. Possible values: {true, false}
maxPopulationAf                      Optional<String>             --max-population-af                                     (-max-af)  Maximum population allele frequency in tumor-only mode.  Default value: 0.01.
maxReadsPerAlignmentStart            Optional<Integer>            --max-reads-per-alignment-start                         Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.  Default value: 50.
minBaseQualityScore                  Optional<String>             --min-base-quality-score                                (-mbq:Byte)  Minimum base quality required to consider a base for calling  Default value: 10.
mitochondriaMode                     Optional<Boolean>            --mitochondria-mode                                     Mitochondria mode sets emission and initial LODs to 0. Default value: false. Possible values: {true, false}
nativePairHmmThreads                 Optional<Integer>            --native-pair-hmm-threads                               How many threads should a native pairHMM implementation use  Default value: 4.
nativePairHmmUseDoublePrecision      Optional<Boolean>            --native-pair-hmm-use-double-precision                  use double precision in the native pairHmm. This is slower but matches the java implementation better  Default value: false. Possible values: {true, false}
normalLod                            Optional<Double>             --normal-lod                                            Log 10 odds threshold for calling normal variant non-germline. Default value: 2.2.
encode                               Optional<String>             -encode                                                 This argument may be specified 0 or more times. Default value: null.
panelOfNormals                       Optional<Gzipped<VCF>>       --panel-of-normals                                      (--panel-of-normals)  VCF file of sites observed in normal.  Default value: null.
pcrIndelQual                         Optional<Integer>            --pcr-indel-qual                                        Phred-scaled PCR SNV qual for overlapping fragments Default value: 40.
pcrSnvQual                           Optional<Integer>            --pcr-snv-qual                                          Phred-scaled PCR SNV qual for overlapping fragments Default value: 40.
pedigree                             Optional<String>             --pedigree                                              (-ped) Pedigree file for determining the population founders. Default value: null.
quiet                                Optional<Boolean>            --QUIET                                                 Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false}
readFilter                           Optional<String>             --read-filter                                           (-RF) Read filters to be applied before analysis This argument may be specified 0 or more times. Default value: null. Possible Values: {AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter}
readIndex                            Optional<String>             -read-index                                             (--read-index)  Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically.  This argument may be specified 0 or more times. Default value: null.
readValidationStringency             Optional<String>             --read-validation-stringency                            (-VS:ValidationStringency)  Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: SILENT. Possible values: {STRICT, LENIENT, SILENT}
secondsBetweenProgressUpdates        Optional<Double>             -seconds-between-progress-updates                       (--seconds-between-progress-updates)  Output traversal statistics every time this many seconds elapse  Default value: 10.0.
sequenceDictionary                   Optional<String>             -sequence-dictionary                                    (--sequence-dictionary)  Use the given sequence dictionary as the master/canonical sequence dictionary.  Must be a .dict file.  Default value: null.
sitesOnlyVcfOutput                   Optional<Boolean>            --sites-only-vcf-output                                 If true, don't emit genotype fields when writing vcf file output.  Default value: false. Possible values: {true, false}
tmpDir                               Optional<String>             --tmp-dir                                               Temp directory to use. Default value: null.
tumorLodToEmit                       Optional<String>             --tumor-lod-to-emit                                     (-emit-lod)  Log 10 odds threshold to emit variant to VCF.  Default value: 3.0.
tumor                                Optional<String>             -tumor                                                  (--tumor-sample) BAM sample name of tumor. May be URL-encoded as output by GetSampleName with -encode argument.  Default value: null.
jdkDeflater                          Optional<Boolean>            -jdk-deflater                                           (--use-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false}
jdkInflater                          Optional<Boolean>            -jdk-inflater                                           (--use-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false}
verbosity                            Optional<String>             -verbosity                                              (--verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG}
version                              Optional<Boolean>            --version                                               display the version number for this tool Default value: false. Possible values: {true, false}
activeProbabilityThreshold           Optional<Double>             --active-probability-threshold                          Minimum probability for a locus to be considered active.  Default value: 0.002.
adaptivePruningInitialErrorRate      Optional<Double>             --adaptive-pruning-initial-error-rate                   Initial base error rate estimate for adaptive pruning  Default value: 0.001.
allowNonUniqueKmersInRef             Optional<Boolean>            --allow-non-unique-kmers-in-ref                         Allow graphs that have non-unique kmers in the reference  Default value: false. Possible values: {true, false}
assemblyRegionPadding                Optional<Integer>            --assembly-region-padding                               Number of additional bases of context to include around each assembly region  Default value: 100.
bamWriterType                        Optional<String>             --bam-writer-type                                       Which haplotypes should be written to the BAM Default value: CALLED_HAPLOTYPES. Possible values: {ALL_POSSIBLE_HAPLOTYPES, CALLED_HAPLOTYPES}
debugAssembly                        Optional<String>             --debug-assembly                                        (-debug)  Print out verbose debug information about each assembly region  Default value: false. Possible values: {true, false}
disableAdaptivePruning               Optional<Boolean>            --disable-adaptive-pruning                              Disable the adaptive algorithm for pruning paths in the graph  Default value: false. Possible values: {true, false}
disableToolDefaultAnnotations        Optional<Boolean>            -disable-tool-default-annotations                       (--disable-tool-default-annotations)  Disable all tool default annotations  Default value: false. Possible values: {true, false}
disableToolDefaultReadFilters        Optional<Boolean>            -disable-tool-default-read-filters                      (--disable-tool-default-read-filters)  Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)  Default value: false. Possible values: {true, false}
dontIncreaseKmerSizesForCycles       Optional<Boolean>            --dont-increase-kmer-sizes-for-cycles                   Disable iterating over kmer sizes when graph cycles are detected  Default value: false. Possible values: {true, false}
dontTrimActiveRegions                Optional<Boolean>            --dont-trim-active-regions                              If specified, we will not trim down the active region from the full region (active + extension) to just the active interval for genotyping  Default value: false. Possible values: {true, false}
dontUseSoftClippedBases              Optional<Boolean>            --dont-use-soft-clipped-bases                           Do not analyze soft clipped bases in the reads  Default value: false. Possible values: {true, false}
erc                                  Optional<String>             -ERC                                                    (--emit-ref-confidence)  (BETA feature) Mode for emitting reference confidence scores  Default value: NONE. Possible values: {NONE, BP_RESOLUTION, GVCF}
enableAllAnnotations                 Optional<Boolean>            --enable-all-annotations                                Use all possible annotations (not for the faint of heart)  Default value: false. Possible values: {true, false}
forceActive                          Optional<Boolean>            --force-active                                          If provided, all regions will be marked as active Default value: false. Possible values: {true, false}
genotypeFilteredAlleles              Optional<Boolean>            --genotype-filtered-alleles                             Whether to force genotype even filtered alleles  Default value: false. Possible values: {true, false}
gvcfLodBand                          Optional<String>             --gvcf-lod-band                                         (-LODB) Exclusive upper bounds for reference confidence LOD bands (must be specified in increasing order)  This argument may be specified 0 or more times. Default value: [-2.5, -2.0, -1.5,
kmerSize                             Optional<Integer>            --kmer-size                                             Kmer size to use in the read threading assembler This argument may be specified 0 or more times. Default value: [10, 25].
maxAssemblyRegionSize                Optional<Integer>            --max-assembly-region-size                              Maximum size of an assembly region  Default value: 300.
mnpDist                              Optional<Integer>            -mnp-dist                                               (--max-mnp-distance)  Two or more phased substitutions separated by this distance or less are merged into MNPs.  Default value: 1.
maxNumHaplotypesInPopulation         Optional<Integer>            --max-num-haplotypes-in-population                      Maximum number of haplotypes to consider for your population  Default value: 128.
maxProbPropagationDistance           Optional<Integer>            --max-prob-propagation-distance                         Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions  Default value: 50.
maxSuspiciousReadsPerAlignmentStart  Optional<Integer>            --max-suspicious-reads-per-alignment-start              Maximum number of suspicious reads (mediocre mapping quality or too many substitutions) allowed in a downsampling stride.  Set to 0 to disable.  Default value: 0.
maxUnprunedVariants                  Optional<Integer>            --max-unpruned-variants                                 Maximum number of variants in graph the adaptive pruner will allow  Default value: 100.
minAssemblyRegionSize                Optional<Integer>            --min-assembly-region-size                              Minimum size of an assembly region  Default value: 50.
minDanglingBranchLength              Optional<Integer>            --min-dangling-branch-length                            Minimum length of a dangling branch to attempt recovery  Default value: 4.
minPruning                           Optional<Integer>            --min-pruning                                           Minimum support to not prune paths in the graph Default value: 2.
minimumAlleleFraction                Optional<Float>              --minimum-allele-fraction                               (-min-AF)  Lower bound of variant allele fractions to consider when calculating variant LOD  Default value: 0.0.
numPruningSamples                    Optional<Integer>            --num-pruning-samples                                   Default value: 1.
pairHmmGapContinuationPenalty        Optional<Integer>            --pair-hmm-gap-continuation-penalty                     Flat gap continuation penalty for use in the Pair HMM  Default value: 10.
pairhmm                              Optional<String>             -pairHMM                                                (--pair-hmm-implementation)  The PairHMM implementation to use for genotype likelihood calculations  Default value: FASTEST_AVAILABLE. Possible values: {EXACT, ORIGINAL, LOGLESS_CACHING, AVX_LOGLESS_CACHING, AVX_LOGLESS_CACHING_OMP, EXPERIMENTAL_FPGA_LOGLESS_CACHING, FASTEST_AVAILABLE}
pcrIndelModel                        Optional<String>             --pcr-indel-model                                       The PCR indel model to use  Default value: CONSERVATIVE. Possible values: {NONE, HOSTILE, AGGRESSIVE, CONSERVATIVE}
phredScaledGlobalReadMismappingRate  Optional<Integer>            --phred-scaled-global-read-mismapping-rate              The global assumed mismapping rate for reads  Default value: 45.
pruningLodThreshold                  Optional<Float>              --pruning-lod-thresholdLn                               Default value: 2.302585092994046.
recoverAllDanglingBranches           Optional<Boolean>            --recover-all-dangling-branches                         Recover all dangling branches  Default value: false. Possible values: {true, false}
showhidden                           Optional<Boolean>            -showHidden                                             (--showHidden)  display hidden arguments  Default value: false. Possible values: {true, false}
smithWaterman                        Optional<String>             --smith-waterman                                        Which Smith-Waterman implementation to use, generally FASTEST_AVAILABLE is the right choice  Default value: JAVA. Possible values: {FASTEST_AVAILABLE, AVX_ENABLED, JAVA}
ambigFilterBases                     Optional<Integer>            --ambig-filter-bases                                    Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction.  Default value: null.  Cannot be used in conjuction with argument(s) maxAmbiguousBaseFraction
ambigFilterFrac                      Optional<Double>             --ambig-filter-frac                                     Threshold fraction of ambiguous bases Default value: 0.05. Cannot be used in conjuction with argument(s) maxAmbiguousBases
maxFragmentLength                    Optional<Integer>            --max-fragment-length                                   Default value: 1000000.
minFragmentLength                    Optional<Integer>            --min-fragment-length                                   Default value: 0.
keepIntervals                        Optional<String>             --keep-intervals                                        One or more genomic intervals to keep This argument must be specified at least once. Required.
library                              Optional<String>             -library                                                (--library) Name of the library to keep This argument must be specified at least once. Required.
maximumMappingQuality                Optional<Integer>            --maximum-mapping-quality                               Maximum mapping quality to keep (inclusive)  Default value: null.
minimumMappingQuality                Optional<Integer>            --minimum-mapping-quality                               Minimum mapping quality to keep (inclusive)  Default value: 20.
dontRequireSoftClipsBothEnds         Optional<Boolean>            --dont-require-soft-clips-both-ends                     Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block  Default value: false. Possible values: {true, false}
filterTooShort                       Optional<Integer>            --filter-too-short                                      Minimum number of aligned bases Default value: 30.
platformFilterName                   Optional<String>             --platform-filter-name                                  This argument must be specified at least once. Required.
blackListedLanes                     Optional<String>             --black-listed-lanes                                    Platform unit (PU) to filter out This argument must be specified at least once. Required.
readGroupBlackList                   Optional<String>             --read-group-black-listThe                              This argument must be specified at least once. Required.
keepReadGroup                        Optional<String>             --keep-read-group                                       The name of the read group to keep Required.
maxReadLength                        Optional<Integer>            --max-read-length                                       Keep only reads with length at most equal to the specified value Default value: 2147483647.
minReadLength                        Optional<Integer>            --min-read-length                                       Keep only reads with length at least equal to the specified value Default value: 30.
readName                             Optional<String>             --read-name                                             Keep only reads with this read name Required.
keepReverseStrandOnly                Optional<Boolean>            --keep-reverse-strand-only                              Keep only reads on the reverse strand  Required. Possible values: {true, false}
sample                               Optional<String>             -sample                                                 (--sample) The name of the sample(s) to keep, filtering out all others This argument must be specified at least once. Required.
===================================  ===========================  ==========================================  ==========  ========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4Mutect2 {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       Array[File] tumorBams
       Array[File] tumorBams_bai
       Array[File]? normalBams
       Array[File]? normalBams_bai
       String? normalSample
       String? outputPrefix
       String? outputFilename
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       String? outputBamName
       String? activityProfileOut
       Boolean? addOutputSamProgramRecord
       Boolean? addOutputVcfCommandLine
       String? afOfAllelesNotInResource
       String? alleles
       String? annotation
       String? annotationGroup
       String? annotationsToExclude
       File? arguments_file
       String? assemblyRegionOut
       Int? baseQualityScoreThreshold
       Int? callableDepth
       Int? cloudIndexPrefetchBuffer
       Int? cloudPrefetchBuffer
       Boolean? createOutputBamIndex
       Boolean? createOutputBamMd5
       Boolean? createOutputVariantIndex
       Boolean? createOutputVariantMd5
       Boolean? disableBamIndexCaching
       Boolean? disableReadFilter
       Boolean? disableSequenceDictionaryValidation
       Int? downsamplingStride
       Boolean? excludeIntervals
       Int? f1r2MaxDepth
       Int? f1r2MedianMq
       Int? f1r2MinBq
       String? f1r2TarGz_outputFilename
       String? founderId
       String? gatkConfigFile
       Int? gcsRetries
       String? gcsProjectForRequesterPays
       Boolean? genotypeGermlineSites
       Boolean? genotypePonSites
       File? germlineResource
       File? germlineResource_tbi
       String? graph
       Boolean? help
       String? ignoreItrArtifacts
       String? initialTumorLod
       String? intervalExclusionPadding
       String? imr
       String? ip
       String? isr
       File? intervals
       Boolean? le
       String? maxPopulationAf
       Int? maxReadsPerAlignmentStart
       String? minBaseQualityScore
       Boolean? mitochondriaMode
       Int? nativePairHmmThreads
       Boolean? nativePairHmmUseDoublePrecision
       Float? normalLod
       String? encode
       File? panelOfNormals
       File? panelOfNormals_tbi
       Int? pcrIndelQual
       Int? pcrSnvQual
       String? pedigree
       Boolean? quiet
       String? readFilter
       String? readIndex
       String? readValidationStringency
       Float? secondsBetweenProgressUpdates
       String? sequenceDictionary
       Boolean? sitesOnlyVcfOutput
       String? tmpDir
       String? tumorLodToEmit
       String? tumor
       Boolean? jdkDeflater
       Boolean? jdkInflater
       String? verbosity
       Boolean? version
       Float? activeProbabilityThreshold
       Float? adaptivePruningInitialErrorRate
       Boolean? allowNonUniqueKmersInRef
       Int? assemblyRegionPadding
       String? bamWriterType
       String? debugAssembly
       Boolean? disableAdaptivePruning
       Boolean? disableToolDefaultAnnotations
       Boolean? disableToolDefaultReadFilters
       Boolean? dontIncreaseKmerSizesForCycles
       Boolean? dontTrimActiveRegions
       Boolean? dontUseSoftClippedBases
       String? erc
       Boolean? enableAllAnnotations
       Boolean? forceActive
       Boolean? genotypeFilteredAlleles
       String? gvcfLodBand
       Int? kmerSize
       Int? maxAssemblyRegionSize
       Int? mnpDist
       Int? maxNumHaplotypesInPopulation
       Int? maxProbPropagationDistance
       Int? maxSuspiciousReadsPerAlignmentStart
       Int? maxUnprunedVariants
       Int? minAssemblyRegionSize
       Int? minDanglingBranchLength
       Int? minPruning
       Float? minimumAlleleFraction
       Int? numPruningSamples
       Int? pairHmmGapContinuationPenalty
       String? pairhmm
       String? pcrIndelModel
       Int? phredScaledGlobalReadMismappingRate
       Float? pruningLodThreshold
       Boolean? recoverAllDanglingBranches
       Boolean? showhidden
       String? smithWaterman
       Int? ambigFilterBases
       Float? ambigFilterFrac
       Int? maxFragmentLength
       Int? minFragmentLength
       String? keepIntervals
       String? library
       Int? maximumMappingQuality
       Int? minimumMappingQuality
       Boolean? dontRequireSoftClipsBothEnds
       Int? filterTooShort
       String? platformFilterName
       String? blackListedLanes
       String? readGroupBlackList
       String? keepReadGroup
       Int? maxReadLength
       Int? minReadLength
       String? readName
       Boolean? keepReverseStrandOnly
       String? sample
     }
     command <<<
       set -e
       gatk Mutect2 \
         --java-options '-Xmx~{((select_first([runtime_memory, 16, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         ~{if length(tumorBams) > 0 then "-I '" + sep("' -I '", tumorBams) + "'" else ""} \
         ~{if (defined(normalBams) && length(select_first([normalBams])) > 0) then "-I '" + sep("' -I '", select_first([normalBams])) + "'" else ""} \
         ~{if defined(normalSample) then ("--normal-sample '" + normalSample + "'") else ""} \
         --reference '~{reference}' \
         ~{if defined(outputBamName) then ("-bamout '" + outputBamName + "'") else ""} \
         ~{if defined(activityProfileOut) then ("--activity-profile-out '" + activityProfileOut + "'") else ""} \
         ~{if (defined(addOutputSamProgramRecord) && select_first([addOutputSamProgramRecord])) then "-add-output-sam-program-record" else ""} \
         ~{if (defined(addOutputVcfCommandLine) && select_first([addOutputVcfCommandLine])) then "-add-output-vcf-command-line" else ""} \
         ~{if defined(afOfAllelesNotInResource) then ("--af-of-alleles-not-in-resource '" + afOfAllelesNotInResource + "'") else ""} \
         ~{if defined(alleles) then ("--alleles '" + alleles + "'") else ""} \
         ~{if defined(annotation) then ("--annotation '" + annotation + "'") else ""} \
         ~{if defined(annotationGroup) then ("--annotation-group '" + annotationGroup + "'") else ""} \
         ~{if defined(annotationsToExclude) then ("--annotations-to-exclude '" + annotationsToExclude + "'") else ""} \
         ~{if defined(arguments_file) then ("--arguments_file '" + arguments_file + "'") else ""} \
         ~{if defined(assemblyRegionOut) then ("--assembly-region-out '" + assemblyRegionOut + "'") else ""} \
         ~{if defined(baseQualityScoreThreshold) then ("--base-quality-score-threshold " + baseQualityScoreThreshold) else ''} \
         ~{if defined(callableDepth) then ("--callable-depth " + callableDepth) else ''} \
         ~{if defined(cloudIndexPrefetchBuffer) then ("--cloud-index-prefetch-buffer " + cloudIndexPrefetchBuffer) else ''} \
         ~{if defined(cloudPrefetchBuffer) then ("--cloud-prefetch-buffer " + cloudPrefetchBuffer) else ''} \
         ~{if (defined(createOutputBamIndex) && select_first([createOutputBamIndex])) then "--create-output-bam-index" else ""} \
         ~{if (defined(createOutputBamMd5) && select_first([createOutputBamMd5])) then "--create-output-bam-md5" else ""} \
         ~{if (defined(createOutputVariantIndex) && select_first([createOutputVariantIndex])) then "--create-output-variant-index" else ""} \
         ~{if (defined(createOutputVariantMd5) && select_first([createOutputVariantMd5])) then "--create-output-variant-md5" else ""} \
         ~{if (defined(disableBamIndexCaching) && select_first([disableBamIndexCaching])) then "--disable-bam-index-caching" else ""} \
         ~{if (defined(disableReadFilter) && select_first([disableReadFilter])) then "--disable-read-filter" else ""} \
         ~{if (defined(disableSequenceDictionaryValidation) && select_first([disableSequenceDictionaryValidation])) then "-disable-sequence-dictionary-validation" else ""} \
         ~{if defined(downsamplingStride) then ("--downsampling-stride " + downsamplingStride) else ''} \
         ~{if (defined(excludeIntervals) && select_first([excludeIntervals])) then "--exclude-intervals" else ""} \
         ~{if defined(f1r2MaxDepth) then ("--f1r2-max-depth " + f1r2MaxDepth) else ''} \
         ~{if defined(f1r2MedianMq) then ("--f1r2-median-mq " + f1r2MedianMq) else ''} \
         ~{if defined(f1r2MinBq) then ("--f1r2-min-bq " + f1r2MinBq) else ''} \
         --f1r2-tar-gz '~{select_first([f1r2TarGz_outputFilename, "generated.tar.gz"])}' \
         ~{if defined(founderId) then ("-founder-id '" + founderId + "'") else ""} \
         ~{if defined(gatkConfigFile) then ("--gatk-config-file '" + gatkConfigFile + "'") else ""} \
         ~{if defined(gcsRetries) then ("-gcs-retries " + gcsRetries) else ''} \
         ~{if defined(gcsProjectForRequesterPays) then ("--gcs-project-for-requester-pays '" + gcsProjectForRequesterPays + "'") else ""} \
         ~{if (defined(genotypeGermlineSites) && select_first([genotypeGermlineSites])) then "--genotype-germline-sites" else ""} \
         ~{if (defined(genotypePonSites) && select_first([genotypePonSites])) then "--genotype-pon-sites" else ""} \
         ~{if defined(germlineResource) then ("--germline-resource '" + germlineResource + "'") else ""} \
         ~{if defined(graph) then ("-graph '" + graph + "'") else ""} \
         ~{if (defined(help) && select_first([help])) then "-h" else ""} \
         ~{if defined(ignoreItrArtifacts) then ("--ignore-itr-artifactsTurn '" + ignoreItrArtifacts + "'") else ""} \
         ~{if defined(initialTumorLod) then ("--initial-tumor-lod '" + initialTumorLod + "'") else ""} \
         ~{if defined(intervalExclusionPadding) then ("--interval-exclusion-padding '" + intervalExclusionPadding + "'") else ""} \
         ~{if defined(imr) then ("--interval-merging-rule '" + imr + "'") else ""} \
         ~{if defined(ip) then ("-ipAmount '" + ip + "'") else ""} \
         ~{if defined(isr) then ("--interval-set-rule '" + isr + "'") else ""} \
         ~{if defined(intervals) then ("--intervals '" + intervals + "'") else ""} \
         ~{if (defined(le) && select_first([le])) then "-LE" else ""} \
         ~{if defined(maxPopulationAf) then ("--max-population-af '" + maxPopulationAf + "'") else ""} \
         ~{if defined(maxReadsPerAlignmentStart) then ("--max-reads-per-alignment-start " + maxReadsPerAlignmentStart) else ''} \
         ~{if defined(minBaseQualityScore) then ("--min-base-quality-score '" + minBaseQualityScore + "'") else ""} \
         ~{if (defined(mitochondriaMode) && select_first([mitochondriaMode])) then "--mitochondria-mode" else ""} \
         ~{if defined(select_first([nativePairHmmThreads, select_first([runtime_cpu, 1])])) then ("--native-pair-hmm-threads " + select_first([nativePairHmmThreads, select_first([runtime_cpu, 1])])) else ''} \
         ~{if (defined(nativePairHmmUseDoublePrecision) && select_first([nativePairHmmUseDoublePrecision])) then "--native-pair-hmm-use-double-precision" else ""} \
         ~{if defined(normalLod) then ("--normal-lod " + normalLod) else ''} \
         ~{if defined(encode) then ("-encode '" + encode + "'") else ""} \
         ~{if defined(panelOfNormals) then ("--panel-of-normals '" + panelOfNormals + "'") else ""} \
         ~{if defined(pcrIndelQual) then ("--pcr-indel-qual " + pcrIndelQual) else ''} \
         ~{if defined(pcrSnvQual) then ("--pcr-snv-qual " + pcrSnvQual) else ''} \
         ~{if defined(pedigree) then ("--pedigree '" + pedigree + "'") else ""} \
         ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
         ~{if defined(readFilter) then ("--read-filter '" + readFilter + "'") else ""} \
         ~{if defined(readIndex) then ("-read-index '" + readIndex + "'") else ""} \
         ~{if defined(readValidationStringency) then ("--read-validation-stringency '" + readValidationStringency + "'") else ""} \
         ~{if defined(secondsBetweenProgressUpdates) then ("-seconds-between-progress-updates " + secondsBetweenProgressUpdates) else ''} \
         ~{if defined(sequenceDictionary) then ("-sequence-dictionary '" + sequenceDictionary + "'") else ""} \
         ~{if (defined(sitesOnlyVcfOutput) && select_first([sitesOnlyVcfOutput])) then "--sites-only-vcf-output" else ""} \
         ~{if defined(tmpDir) then ("--tmp-dir '" + tmpDir + "'") else ""} \
         ~{if defined(tumorLodToEmit) then ("--tumor-lod-to-emit '" + tumorLodToEmit + "'") else ""} \
         ~{if defined(tumor) then ("-tumor '" + tumor + "'") else ""} \
         ~{if (defined(jdkDeflater) && select_first([jdkDeflater])) then "-jdk-deflater" else ""} \
         ~{if (defined(jdkInflater) && select_first([jdkInflater])) then "-jdk-inflater" else ""} \
         ~{if defined(verbosity) then ("-verbosity '" + verbosity + "'") else ""} \
         ~{if (defined(version) && select_first([version])) then "--version" else ""} \
         ~{if defined(activeProbabilityThreshold) then ("--active-probability-threshold " + activeProbabilityThreshold) else ''} \
         ~{if defined(adaptivePruningInitialErrorRate) then ("--adaptive-pruning-initial-error-rate " + adaptivePruningInitialErrorRate) else ''} \
         ~{if (defined(allowNonUniqueKmersInRef) && select_first([allowNonUniqueKmersInRef])) then "--allow-non-unique-kmers-in-ref" else ""} \
         ~{if defined(assemblyRegionPadding) then ("--assembly-region-padding " + assemblyRegionPadding) else ''} \
         ~{if defined(bamWriterType) then ("--bam-writer-type '" + bamWriterType + "'") else ""} \
         ~{if defined(debugAssembly) then ("--debug-assembly '" + debugAssembly + "'") else ""} \
         ~{if (defined(disableAdaptivePruning) && select_first([disableAdaptivePruning])) then "--disable-adaptive-pruning" else ""} \
         ~{if (defined(disableToolDefaultAnnotations) && select_first([disableToolDefaultAnnotations])) then "-disable-tool-default-annotations" else ""} \
         ~{if (defined(disableToolDefaultReadFilters) && select_first([disableToolDefaultReadFilters])) then "-disable-tool-default-read-filters" else ""} \
         ~{if (defined(dontIncreaseKmerSizesForCycles) && select_first([dontIncreaseKmerSizesForCycles])) then "--dont-increase-kmer-sizes-for-cycles" else ""} \
         ~{if (defined(dontTrimActiveRegions) && select_first([dontTrimActiveRegions])) then "--dont-trim-active-regions" else ""} \
         ~{if (defined(dontUseSoftClippedBases) && select_first([dontUseSoftClippedBases])) then "--dont-use-soft-clipped-bases" else ""} \
         ~{if defined(erc) then ("-ERC '" + erc + "'") else ""} \
         ~{if (defined(enableAllAnnotations) && select_first([enableAllAnnotations])) then "--enable-all-annotations" else ""} \
         ~{if (defined(forceActive) && select_first([forceActive])) then "--force-active" else ""} \
         ~{if (defined(genotypeFilteredAlleles) && select_first([genotypeFilteredAlleles])) then "--genotype-filtered-alleles" else ""} \
         ~{if defined(gvcfLodBand) then ("--gvcf-lod-band '" + gvcfLodBand + "'") else ""} \
         ~{if defined(kmerSize) then ("--kmer-size " + kmerSize) else ''} \
         ~{if defined(maxAssemblyRegionSize) then ("--max-assembly-region-size " + maxAssemblyRegionSize) else ''} \
         ~{if defined(mnpDist) then ("-mnp-dist " + mnpDist) else ''} \
         ~{if defined(maxNumHaplotypesInPopulation) then ("--max-num-haplotypes-in-population " + maxNumHaplotypesInPopulation) else ''} \
         ~{if defined(maxProbPropagationDistance) then ("--max-prob-propagation-distance " + maxProbPropagationDistance) else ''} \
         ~{if defined(maxSuspiciousReadsPerAlignmentStart) then ("--max-suspicious-reads-per-alignment-start " + maxSuspiciousReadsPerAlignmentStart) else ''} \
         ~{if defined(maxUnprunedVariants) then ("--max-unpruned-variants " + maxUnprunedVariants) else ''} \
         ~{if defined(minAssemblyRegionSize) then ("--min-assembly-region-size " + minAssemblyRegionSize) else ''} \
         ~{if defined(minDanglingBranchLength) then ("--min-dangling-branch-length " + minDanglingBranchLength) else ''} \
         ~{if defined(minPruning) then ("--min-pruning " + minPruning) else ''} \
         ~{if defined(minimumAlleleFraction) then ("--minimum-allele-fraction " + minimumAlleleFraction) else ''} \
         ~{if defined(numPruningSamples) then ("--num-pruning-samples " + numPruningSamples) else ''} \
         ~{if defined(pairHmmGapContinuationPenalty) then ("--pair-hmm-gap-continuation-penalty " + pairHmmGapContinuationPenalty) else ''} \
         ~{if defined(pairhmm) then ("-pairHMM '" + pairhmm + "'") else ""} \
         ~{if defined(pcrIndelModel) then ("--pcr-indel-model '" + pcrIndelModel + "'") else ""} \
         ~{if defined(phredScaledGlobalReadMismappingRate) then ("--phred-scaled-global-read-mismapping-rate " + phredScaledGlobalReadMismappingRate) else ''} \
         ~{if defined(pruningLodThreshold) then ("--pruning-lod-thresholdLn " + pruningLodThreshold) else ''} \
         ~{if (defined(recoverAllDanglingBranches) && select_first([recoverAllDanglingBranches])) then "--recover-all-dangling-branches" else ""} \
         ~{if (defined(showhidden) && select_first([showhidden])) then "-showHidden" else ""} \
         ~{if defined(smithWaterman) then ("--smith-waterman '" + smithWaterman + "'") else ""} \
         ~{if defined(ambigFilterBases) then ("--ambig-filter-bases " + ambigFilterBases) else ''} \
         ~{if defined(ambigFilterFrac) then ("--ambig-filter-frac " + ambigFilterFrac) else ''} \
         ~{if defined(maxFragmentLength) then ("--max-fragment-length " + maxFragmentLength) else ''} \
         ~{if defined(minFragmentLength) then ("--min-fragment-length " + minFragmentLength) else ''} \
         ~{if defined(keepIntervals) then ("--keep-intervals '" + keepIntervals + "'") else ""} \
         ~{if defined(library) then ("-library '" + library + "'") else ""} \
         ~{if defined(maximumMappingQuality) then ("--maximum-mapping-quality " + maximumMappingQuality) else ''} \
         ~{if defined(minimumMappingQuality) then ("--minimum-mapping-quality " + minimumMappingQuality) else ''} \
         ~{if (defined(dontRequireSoftClipsBothEnds) && select_first([dontRequireSoftClipsBothEnds])) then "--dont-require-soft-clips-both-ends" else ""} \
         ~{if defined(filterTooShort) then ("--filter-too-short " + filterTooShort) else ''} \
         ~{if defined(platformFilterName) then ("--platform-filter-name '" + platformFilterName + "'") else ""} \
         ~{if defined(blackListedLanes) then ("--black-listed-lanes '" + blackListedLanes + "'") else ""} \
         ~{if defined(readGroupBlackList) then ("--read-group-black-listThe '" + readGroupBlackList + "'") else ""} \
         ~{if defined(keepReadGroup) then ("--keep-read-group '" + keepReadGroup + "'") else ""} \
         ~{if defined(maxReadLength) then ("--max-read-length " + maxReadLength) else ''} \
         ~{if defined(minReadLength) then ("--min-read-length " + minReadLength) else ''} \
         ~{if defined(readName) then ("--read-name '" + readName + "'") else ""} \
         ~{if (defined(keepReverseStrandOnly) && select_first([keepReverseStrandOnly])) then "--keep-reverse-strand-only" else ""} \
         ~{if defined(sample) then ("-sample '" + sample + "'") else ""} \
         -O '~{select_first([outputFilename, "~{select_first([outputPrefix, "generated"])}.vcf.gz"])}'
       if [ -f $(echo '~{outputBamName}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{outputBamName}' | sed 's/\.[^.]*$//').bai $(echo '~{outputBamName}' ).bai; fi
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.4.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 16, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "~{select_first([outputPrefix, "generated"])}.vcf.gz"])
       File out_tbi = select_first([outputFilename, "~{select_first([outputPrefix, "generated"])}.vcf.gz"]) + ".tbi"
       File stats = (select_first([outputFilename, "~{select_first([outputPrefix, "generated"])}.vcf.gz"]) + ".stats")
       File f1f2r_out = select_first([f1r2TarGz_outputFilename, "generated.tar.gz"])
       File? bam = outputBamName
       File? bam_bai = if defined(outputBamName) then (outputBamName + ".bai") else None
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: GatkMutect2
   doc: |
     USAGE: Mutect2 [arguments]
     Call somatic SNVs and indels via local assembly of haplotypes
     Version:4.1.2.0

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.4.0

   inputs:
   - id: javaOptions
     label: javaOptions
     type:
     - type: array
       items: string
     - 'null'
   - id: compression_level
     label: compression_level
     doc: |-
       Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
     type:
     - int
     - 'null'
   - id: tumorBams
     label: tumorBams
     doc: |-
       (--input) BAM/SAM/CRAM file containing reads This argument must be specified at least once. Required. 
     type:
       type: array
       inputBinding:
         prefix: -I
       items: File
     inputBinding: {}
   - id: normalBams
     label: normalBams
     doc: |-
       (--input) Extra BAM/SAM/CRAM file containing reads This argument must be specified at least once. Required. 
     type:
     - type: array
       inputBinding:
         prefix: -I
       items: File
     - 'null'
     inputBinding: {}
   - id: normalSample
     label: normalSample
     doc: (--normal-sample, if) May be URL-encoded as output by GetSampleName with
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --normal-sample
   - id: outputPrefix
     label: outputPrefix
     doc: |-
       Used as a prefix for the outputFilename if not specified, with format: {outputPrefix}.vcf.gz
     type: string
     default: generated
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.vcf.gz
     inputBinding:
       prefix: -O
       position: 20
       valueFrom: '$(inputs.outputPrefix ? inputs.outputPrefix : "generated").vcf.gz'
   - id: reference
     label: reference
     doc: (-R) Reference sequence file Required.
     type: File
     secondaryFiles:
     - pattern: .fai
     - pattern: .amb
     - pattern: .ann
     - pattern: .bwt
     - pattern: .pac
     - pattern: .sa
     - pattern: ^.dict
     inputBinding:
       prefix: --reference
   - id: outputBamName
     label: outputBamName
     doc: File to which assembled haplotypes should be written
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -bamout
   - id: activityProfileOut
     label: activityProfileOut
     doc: 'Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --activity-profile-out
   - id: addOutputSamProgramRecord
     label: addOutputSamProgramRecord
     doc: |-
       (--add-output-sam-program-record)  If true, adds a PG tag to created SAM/BAM/CRAM files.  Default value: true. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -add-output-sam-program-record
   - id: addOutputVcfCommandLine
     label: addOutputVcfCommandLine
     doc: |-
       (--add-output-vcf-command-line)  If true, adds a command line header line to created VCF files.  Default value: true. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -add-output-vcf-command-line
   - id: afOfAllelesNotInResource
     label: afOfAllelesNotInResource
     doc: |-
       (-default-af)  Population allele fraction assigned to alleles not found in germline resource.  Please see docs/mutect/mutect2.pdf fora derivation of the default value.  Default value: -1.0. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --af-of-alleles-not-in-resource
   - id: alleles
     label: alleles
     doc: |-
       The set of alleles for which to force genotyping regardless of evidence Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --alleles
   - id: annotation
     label: annotation
     doc: |-
       (-A) One or more specific annotations to add to variant calls This argument may be specified 0 or more times. Default value: null. Possible Values: {AlleleFraction, AS_BaseQualityRankSumTest, AS_FisherStrand, AS_InbreedingCoeff, AS_MappingQualityRankSumTest, AS_QualByDepth, AS_ReadPosRankSumTest, AS_RMSMappingQuality, AS_StrandOddsRatio, BaseQuality, BaseQualityRankSumTest, ChromosomeCounts, ClippingRankSumTest, CountNs, Coverage, DepthPerAlleleBySample, DepthPerSampleHC, ExcessHet, FisherStrand, FragmentLength, GenotypeSummaries, InbreedingCoeff, LikelihoodRankSumTest, MappingQuality, MappingQualityRankSumTest, MappingQualityZero, OrientationBiasReadCounts, OriginalAlignment, PossibleDeNovo, QualByDepth, ReadPosition, ReadPosRankSumTest, ReferenceBases, RMSMappingQuality, SampleList, StrandBiasBySample, StrandOddsRatio, TandemRepeat, UniqueAltReadCount}
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --annotation
   - id: annotationGroup
     label: annotationGroup
     doc: |-
       (-G) One or more groups of annotations to apply to variant calls This argument may be specified 0 or more times. Default value: null. Possible Values: {AS_StandardAnnotation, ReducibleAnnotation, StandardAnnotation, StandardHCAnnotation, StandardMutectAnnotation}
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --annotation-group
   - id: annotationsToExclude
     label: annotationsToExclude
     doc: |-
       (-AX)  One or more specific annotations to exclude from variant calls  This argument may be specified 0 or more times. Default value: null. Possible Values: {BaseQuality, Coverage, DepthPerAlleleBySample, DepthPerSampleHC, FragmentLength, MappingQuality, OrientationBiasReadCounts, ReadPosition, StrandBiasBySample, TandemRepeat}
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --annotations-to-exclude
   - id: arguments_file
     label: arguments_file
     doc: |-
       read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null. 
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --arguments_file
   - id: assemblyRegionOut
     label: assemblyRegionOut
     doc: 'Output the assembly region to this IGV formatted file Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --assembly-region-out
   - id: baseQualityScoreThreshold
     label: baseQualityScoreThreshold
     doc: |2-
        Base qualities below this threshold will be reduced to the minimum (6)  Default value: 18.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --base-quality-score-threshold
   - id: callableDepth
     label: callableDepth
     doc: |-
       Minimum depth to be considered callable for Mutect stats. Does not affect genotyping. Default value: 10. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --callable-depth
   - id: cloudIndexPrefetchBuffer
     label: cloudIndexPrefetchBuffer
     doc: |-
       (-CIPB)  Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.  Default value: -1. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --cloud-index-prefetch-buffer
   - id: cloudPrefetchBuffer
     label: cloudPrefetchBuffer
     doc: |-
       (-CPB)  Size of the cloud-only prefetch buffer (in MB; 0 to disable).  Default value: 40. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --cloud-prefetch-buffer
   - id: createOutputBamIndex
     label: createOutputBamIndex
     doc: |-
       (-OBI)  If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.  Default value: true. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --create-output-bam-index
   - id: createOutputBamMd5
     label: createOutputBamMd5
     doc: |-
       (-OBM)  If true, create a MD5 digest for any BAM/SAM/CRAM file created  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --create-output-bam-md5
   - id: createOutputVariantIndex
     label: createOutputVariantIndex
     doc: |-
       (-OVI)  If true, create a VCF index when writing a coordinate-sorted VCF file.  Default value: true. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --create-output-variant-index
   - id: createOutputVariantMd5
     label: createOutputVariantMd5
     doc: |-
       (-OVM)  If true, create a a MD5 digest any VCF file created.  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --create-output-variant-md5
   - id: disableBamIndexCaching
     label: disableBamIndexCaching
     doc: |-
       (-DBIC)  If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --disable-bam-index-caching
   - id: disableReadFilter
     label: disableReadFilter
     doc: |-
       (-DF)  Read filters to be disabled before analysis  This argument may be specified 0 or more times. Default value: null. Possible Values: {GoodCigarReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotSecondaryAlignmentReadFilter, PassesVendorQualityCheckReadFilter, ReadLengthReadFilter, WellformedReadFilter}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --disable-read-filter
   - id: disableSequenceDictionaryValidation
     label: disableSequenceDictionaryValidation
     doc: |-
       (--disable-sequence-dictionary-validation)  If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -disable-sequence-dictionary-validation
   - id: downsamplingStride
     label: downsamplingStride
     doc: |-
       (-stride)  Downsample a pool of reads starting within a range of one or more bases.  Default value: 1. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --downsampling-stride
   - id: excludeIntervals
     label: excludeIntervals
     doc: '(-XLOne) This argument may be specified 0 or more times. Default value: null. '
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --exclude-intervals
   - id: f1r2MaxDepth
     label: f1r2MaxDepth
     doc: 'sites with depth higher than this value will be grouped Default value: 200.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --f1r2-max-depth
   - id: f1r2MedianMq
     label: f1r2MedianMq
     doc: 'skip sites with median mapping quality below this value Default value: 50.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --f1r2-median-mq
   - id: f1r2MinBq
     label: f1r2MinBq
     doc: 'exclude bases below this quality from pileup Default value: 20.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --f1r2-min-bq
   - id: f1r2TarGz_outputFilename
     label: f1r2TarGz_outputFilename
     doc: |-
       If specified, collect F1R2 counts and output files into this tar.gz file Default value: null. 
     type:
     - string
     - 'null'
     default: generated.tar.gz
     inputBinding:
       prefix: --f1r2-tar-gz
   - id: founderId
     label: founderId
     doc: |-
       (--founder-id)  Samples representing the population founders This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -founder-id
   - id: gatkConfigFile
     label: gatkConfigFile
     doc: 'A configuration file to use with the GATK. Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --gatk-config-file
   - id: gcsRetries
     label: gcsRetries
     doc: |-
       (--gcs-max-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -gcs-retries
   - id: gcsProjectForRequesterPays
     label: gcsProjectForRequesterPays
     doc: |2-
        Project to bill when accessing requester pays buckets. If unset, these buckets cannot be accessed.  Default value: . 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --gcs-project-for-requester-pays
   - id: genotypeGermlineSites
     label: genotypeGermlineSites
     doc: |2-
        (EXPERIMENTAL) Call all apparent germline site even though they will ultimately be filtered.  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --genotype-germline-sites
   - id: genotypePonSites
     label: genotypePonSites
     doc: |-
       Call sites in the PoN even though they will ultimately be filtered. Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --genotype-pon-sites
   - id: germlineResource
     label: germlineResource
     doc: |2-
        Population vcf of germline sequencing containing allele fractions.  Default value: null. 
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       prefix: --germline-resource
   - id: graph
     label: graph
     doc: |-
       (--graph-output) Write debug assembly graph information to this file Default value: null.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -graph
   - id: help
     label: help
     doc: |-
       (--help) display the help message Default value: false. Possible values: {true, false}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -h
   - id: ignoreItrArtifacts
     label: ignoreItrArtifacts
     doc: ' inverted tandem repeats.  Default value: false. Possible values: {true, false} '
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --ignore-itr-artifactsTurn
   - id: initialTumorLod
     label: initialTumorLod
     doc: |-
       (-init-lod)  Log 10 odds threshold to consider pileup active.  Default value: 2.0. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --initial-tumor-lod
   - id: intervalExclusionPadding
     label: intervalExclusionPadding
     doc: |-
       (-ixp)  Amount of padding (in bp) to add to each interval you are excluding.  Default value: 0. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --interval-exclusion-padding
   - id: imr
     label: imr
     doc: |-
       (--interval-merging-rule)  Interval merging rule for abutting intervals  Default value: ALL. Possible values: {ALL, OVERLAPPING_ONLY} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --interval-merging-rule
   - id: ip
     label: ip
     doc: '(--interval-padding) Default value: 0.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -ipAmount
   - id: isr
     label: isr
     doc: |-
       (--interval-set-rule)  Set merging approach to use for combining interval inputs  Default value: UNION. Possible values: {UNION, INTERSECTION} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --interval-set-rule
   - id: intervals
     label: intervals
     doc: |-
       (-L) One or more genomic intervals over which to operate This argument may be specified 0 or more times. Default value: null. 
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --intervals
   - id: le
     label: le
     doc: |-
       (--lenient) Lenient processing of VCF files Default value: false. Possible values: {true, false}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -LE
   - id: maxPopulationAf
     label: maxPopulationAf
     doc: |-
       (-max-af)  Maximum population allele frequency in tumor-only mode.  Default value: 0.01. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --max-population-af
   - id: maxReadsPerAlignmentStart
     label: maxReadsPerAlignmentStart
     doc: |2-
        Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.  Default value: 50. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-reads-per-alignment-start
   - id: minBaseQualityScore
     label: minBaseQualityScore
     doc: |-
       (-mbq:Byte)  Minimum base quality required to consider a base for calling  Default value: 10. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --min-base-quality-score
   - id: mitochondriaMode
     label: mitochondriaMode
     doc: |-
       Mitochondria mode sets emission and initial LODs to 0. Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --mitochondria-mode
   - id: nativePairHmmThreads
     label: nativePairHmmThreads
     doc: ' How many threads should a native pairHMM implementation use  Default value:
       4. '
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --native-pair-hmm-threads
       valueFrom: |-
         $([inputs.runtime_cpu, 4, 1].filter(function (inner) { return inner != null })[0])
   - id: nativePairHmmUseDoublePrecision
     label: nativePairHmmUseDoublePrecision
     doc: |2-
        use double precision in the native pairHmm. This is slower but matches the java implementation better  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --native-pair-hmm-use-double-precision
   - id: normalLod
     label: normalLod
     doc: |-
       Log 10 odds threshold for calling normal variant non-germline. Default value: 2.2.
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --normal-lod
   - id: encode
     label: encode
     doc: 'This argument may be specified 0 or more times. Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -encode
   - id: panelOfNormals
     label: panelOfNormals
     doc: |-
       (--panel-of-normals)  VCF file of sites observed in normal.  Default value: null. 
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       prefix: --panel-of-normals
   - id: pcrIndelQual
     label: pcrIndelQual
     doc: 'Phred-scaled PCR SNV qual for overlapping fragments Default value: 40.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --pcr-indel-qual
   - id: pcrSnvQual
     label: pcrSnvQual
     doc: 'Phred-scaled PCR SNV qual for overlapping fragments Default value: 40.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --pcr-snv-qual
   - id: pedigree
     label: pedigree
     doc: |-
       (-ped) Pedigree file for determining the population founders. Default value: null.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --pedigree
   - id: quiet
     label: quiet
     doc: |-
       Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --QUIET
   - id: readFilter
     label: readFilter
     doc: |-
       (-RF) Read filters to be applied before analysis This argument may be specified 0 or more times. Default value: null. Possible Values: {AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter}
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read-filter
   - id: readIndex
     label: readIndex
     doc: |-
       (--read-index)  Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically.  This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -read-index
   - id: readValidationStringency
     label: readValidationStringency
     doc: |-
       (-VS:ValidationStringency)  Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: SILENT. Possible values: {STRICT, LENIENT, SILENT} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read-validation-stringency
   - id: secondsBetweenProgressUpdates
     label: secondsBetweenProgressUpdates
     doc: |-
       (--seconds-between-progress-updates)  Output traversal statistics every time this many seconds elapse  Default value: 10.0. 
     type:
     - double
     - 'null'
     inputBinding:
       prefix: -seconds-between-progress-updates
   - id: sequenceDictionary
     label: sequenceDictionary
     doc: |-
       (--sequence-dictionary)  Use the given sequence dictionary as the master/canonical sequence dictionary.  Must be a .dict file.  Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -sequence-dictionary
   - id: sitesOnlyVcfOutput
     label: sitesOnlyVcfOutput
     doc: |2-
        If true, don't emit genotype fields when writing vcf file output.  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --sites-only-vcf-output
   - id: tmpDir
     label: tmpDir
     doc: 'Temp directory to use. Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --tmp-dir
   - id: tumorLodToEmit
     label: tumorLodToEmit
     doc: '(-emit-lod)  Log 10 odds threshold to emit variant to VCF.  Default value:
       3.0. '
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --tumor-lod-to-emit
   - id: tumor
     label: tumor
     doc: |-
       (--tumor-sample) BAM sample name of tumor. May be URL-encoded as output by GetSampleName with -encode argument.  Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -tumor
   - id: jdkDeflater
     label: jdkDeflater
     doc: |-
       (--use-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -jdk-deflater
   - id: jdkInflater
     label: jdkInflater
     doc: |-
       (--use-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -jdk-inflater
   - id: verbosity
     label: verbosity
     doc: |-
       (--verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -verbosity
   - id: version
     label: version
     doc: |-
       display the version number for this tool Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --version
   - id: activeProbabilityThreshold
     label: activeProbabilityThreshold
     doc: |2-
        Minimum probability for a locus to be considered active.  Default value: 0.002. 
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --active-probability-threshold
   - id: adaptivePruningInitialErrorRate
     label: adaptivePruningInitialErrorRate
     doc: ' Initial base error rate estimate for adaptive pruning  Default value: 0.001. '
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --adaptive-pruning-initial-error-rate
   - id: allowNonUniqueKmersInRef
     label: allowNonUniqueKmersInRef
     doc: |2-
        Allow graphs that have non-unique kmers in the reference  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --allow-non-unique-kmers-in-ref
   - id: assemblyRegionPadding
     label: assemblyRegionPadding
     doc: |2-
        Number of additional bases of context to include around each assembly region  Default value: 100. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --assembly-region-padding
   - id: bamWriterType
     label: bamWriterType
     doc: |-
       Which haplotypes should be written to the BAM Default value: CALLED_HAPLOTYPES. Possible values: {ALL_POSSIBLE_HAPLOTYPES, CALLED_HAPLOTYPES} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --bam-writer-type
   - id: debugAssembly
     label: debugAssembly
     doc: |-
       (-debug)  Print out verbose debug information about each assembly region  Default value: false. Possible values: {true, false} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --debug-assembly
   - id: disableAdaptivePruning
     label: disableAdaptivePruning
     doc: |2-
        Disable the adaptive algorithm for pruning paths in the graph  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --disable-adaptive-pruning
   - id: disableToolDefaultAnnotations
     label: disableToolDefaultAnnotations
     doc: |-
       (--disable-tool-default-annotations)  Disable all tool default annotations  Default value: false. Possible values: {true, false}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -disable-tool-default-annotations
   - id: disableToolDefaultReadFilters
     label: disableToolDefaultReadFilters
     doc: |-
       (--disable-tool-default-read-filters)  Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -disable-tool-default-read-filters
   - id: dontIncreaseKmerSizesForCycles
     label: dontIncreaseKmerSizesForCycles
     doc: |2-
        Disable iterating over kmer sizes when graph cycles are detected  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --dont-increase-kmer-sizes-for-cycles
   - id: dontTrimActiveRegions
     label: dontTrimActiveRegions
     doc: |2-
        If specified, we will not trim down the active region from the full region (active + extension) to just the active interval for genotyping  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --dont-trim-active-regions
   - id: dontUseSoftClippedBases
     label: dontUseSoftClippedBases
     doc: |2-
        Do not analyze soft clipped bases in the reads  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --dont-use-soft-clipped-bases
   - id: erc
     label: erc
     doc: |-
       (--emit-ref-confidence)  (BETA feature) Mode for emitting reference confidence scores  Default value: NONE. Possible values: {NONE, BP_RESOLUTION, GVCF} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -ERC
   - id: enableAllAnnotations
     label: enableAllAnnotations
     doc: |2-
        Use all possible annotations (not for the faint of heart)  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --enable-all-annotations
   - id: forceActive
     label: forceActive
     doc: |-
       If provided, all regions will be marked as active Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --force-active
   - id: genotypeFilteredAlleles
     label: genotypeFilteredAlleles
     doc: |2-
        Whether to force genotype even filtered alleles  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --genotype-filtered-alleles
   - id: gvcfLodBand
     label: gvcfLodBand
     doc: |-
       (-LODB) Exclusive upper bounds for reference confidence LOD bands (must be specified in increasing order)  This argument may be specified 0 or more times. Default value: [-2.5, -2.0, -1.5,
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --gvcf-lod-band
   - id: kmerSize
     label: kmerSize
     doc: |-
       Kmer size to use in the read threading assembler This argument may be specified 0 or more times. Default value: [10, 25]. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --kmer-size
   - id: maxAssemblyRegionSize
     label: maxAssemblyRegionSize
     doc: ' Maximum size of an assembly region  Default value: 300. '
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-assembly-region-size
   - id: mnpDist
     label: mnpDist
     doc: |-
       (--max-mnp-distance)  Two or more phased substitutions separated by this distance or less are merged into MNPs.  Default value: 1. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -mnp-dist
   - id: maxNumHaplotypesInPopulation
     label: maxNumHaplotypesInPopulation
     doc: |2-
        Maximum number of haplotypes to consider for your population  Default value: 128. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-num-haplotypes-in-population
   - id: maxProbPropagationDistance
     label: maxProbPropagationDistance
     doc: |2-
        Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions  Default value: 50. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-prob-propagation-distance
   - id: maxSuspiciousReadsPerAlignmentStart
     label: maxSuspiciousReadsPerAlignmentStart
     doc: |2-
        Maximum number of suspicious reads (mediocre mapping quality or too many substitutions) allowed in a downsampling stride.  Set to 0 to disable.  Default value: 0. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-suspicious-reads-per-alignment-start
   - id: maxUnprunedVariants
     label: maxUnprunedVariants
     doc: |2-
        Maximum number of variants in graph the adaptive pruner will allow  Default value: 100. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-unpruned-variants
   - id: minAssemblyRegionSize
     label: minAssemblyRegionSize
     doc: ' Minimum size of an assembly region  Default value: 50. '
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-assembly-region-size
   - id: minDanglingBranchLength
     label: minDanglingBranchLength
     doc: ' Minimum length of a dangling branch to attempt recovery  Default value: 4. '
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-dangling-branch-length
   - id: minPruning
     label: minPruning
     doc: 'Minimum support to not prune paths in the graph Default value: 2.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-pruning
   - id: minimumAlleleFraction
     label: minimumAlleleFraction
     doc: |-
       (-min-AF)  Lower bound of variant allele fractions to consider when calculating variant LOD  Default value: 0.0. 
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --minimum-allele-fraction
   - id: numPruningSamples
     label: numPruningSamples
     doc: 'Default value: 1.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --num-pruning-samples
   - id: pairHmmGapContinuationPenalty
     label: pairHmmGapContinuationPenalty
     doc: ' Flat gap continuation penalty for use in the Pair HMM  Default value: 10. '
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --pair-hmm-gap-continuation-penalty
   - id: pairhmm
     label: pairhmm
     doc: |-
       (--pair-hmm-implementation)  The PairHMM implementation to use for genotype likelihood calculations  Default value: FASTEST_AVAILABLE. Possible values: {EXACT, ORIGINAL, LOGLESS_CACHING, AVX_LOGLESS_CACHING, AVX_LOGLESS_CACHING_OMP, EXPERIMENTAL_FPGA_LOGLESS_CACHING, FASTEST_AVAILABLE} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -pairHMM
   - id: pcrIndelModel
     label: pcrIndelModel
     doc: |2-
        The PCR indel model to use  Default value: CONSERVATIVE. Possible values: {NONE, HOSTILE, AGGRESSIVE, CONSERVATIVE} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --pcr-indel-model
   - id: phredScaledGlobalReadMismappingRate
     label: phredScaledGlobalReadMismappingRate
     doc: ' The global assumed mismapping rate for reads  Default value: 45. '
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --phred-scaled-global-read-mismapping-rate
   - id: pruningLodThreshold
     label: pruningLodThreshold
     doc: 'Default value: 2.302585092994046. '
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --pruning-lod-thresholdLn
   - id: recoverAllDanglingBranches
     label: recoverAllDanglingBranches
     doc: |2-
        Recover all dangling branches  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --recover-all-dangling-branches
   - id: showhidden
     label: showhidden
     doc: |-
       (--showHidden)  display hidden arguments  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -showHidden
   - id: smithWaterman
     label: smithWaterman
     doc: |2-
        Which Smith-Waterman implementation to use, generally FASTEST_AVAILABLE is the right choice  Default value: JAVA. Possible values: {FASTEST_AVAILABLE, AVX_ENABLED, JAVA} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --smith-waterman
   - id: ambigFilterBases
     label: ambigFilterBases
     doc: |-
       Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction.  Default value: null.  Cannot be used in conjuction with argument(s) maxAmbiguousBaseFraction
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --ambig-filter-bases
   - id: ambigFilterFrac
     label: ambigFilterFrac
     doc: |-
       Threshold fraction of ambiguous bases Default value: 0.05. Cannot be used in conjuction with argument(s) maxAmbiguousBases
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --ambig-filter-frac
   - id: maxFragmentLength
     label: maxFragmentLength
     doc: 'Default value: 1000000.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-fragment-length
   - id: minFragmentLength
     label: minFragmentLength
     doc: 'Default value: 0.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-fragment-length
   - id: keepIntervals
     label: keepIntervals
     doc: |-
       One or more genomic intervals to keep This argument must be specified at least once. Required. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --keep-intervals
   - id: library
     label: library
     doc: |-
       (--library) Name of the library to keep This argument must be specified at least once. Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -library
   - id: maximumMappingQuality
     label: maximumMappingQuality
     doc: ' Maximum mapping quality to keep (inclusive)  Default value: null. '
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --maximum-mapping-quality
   - id: minimumMappingQuality
     label: minimumMappingQuality
     doc: ' Minimum mapping quality to keep (inclusive)  Default value: 20. '
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --minimum-mapping-quality
   - id: dontRequireSoftClipsBothEnds
     label: dontRequireSoftClipsBothEnds
     doc: |2-
        Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --dont-require-soft-clips-both-ends
   - id: filterTooShort
     label: filterTooShort
     doc: 'Minimum number of aligned bases Default value: 30.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --filter-too-short
   - id: platformFilterName
     label: platformFilterName
     doc: This argument must be specified at least once. Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --platform-filter-name
   - id: blackListedLanes
     label: blackListedLanes
     doc: |-
       Platform unit (PU) to filter out This argument must be specified at least once. Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --black-listed-lanes
   - id: readGroupBlackList
     label: readGroupBlackList
     doc: 'This argument must be specified at least once. Required. '
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read-group-black-listThe
   - id: keepReadGroup
     label: keepReadGroup
     doc: The name of the read group to keep Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --keep-read-group
   - id: maxReadLength
     label: maxReadLength
     doc: |-
       Keep only reads with length at most equal to the specified value Default value: 2147483647. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-read-length
   - id: minReadLength
     label: minReadLength
     doc: |-
       Keep only reads with length at least equal to the specified value Default value: 30.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-read-length
   - id: readName
     label: readName
     doc: Keep only reads with this read name Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read-name
   - id: keepReverseStrandOnly
     label: keepReverseStrandOnly
     doc: |2-
        Keep only reads on the reverse strand  Required. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --keep-reverse-strand-only
   - id: sample
     label: sample
     doc: |-
       (--sample) The name of the sample(s) to keep, filtering out all others This argument must be specified at least once. Required. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -sample

   outputs:
   - id: out
     label: out
     doc: To determine type
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: '$(inputs.outputPrefix ? inputs.outputPrefix : "generated").vcf.gz'
       loadContents: false
   - id: stats
     label: stats
     doc: To determine type
     type: File
     outputBinding:
       glob: $((inputs.outputFilename + ".stats"))
       outputEval: $((inputs.outputFilename.basename + ".stats"))
       loadContents: false
   - id: f1f2r_out
     label: f1f2r_out
     doc: To determine type
     type: File
     outputBinding:
       glob: generated.tar.gz
       loadContents: false
   - id: bam
     label: bam
     doc: File to which assembled haplotypes should be written
     type:
     - File
     - 'null'
     secondaryFiles:
     - |-
       ${

               function resolveSecondary(base, secPattern) {
                 if (secPattern[0] == "^") {
                   var spl = base.split(".");
                   var endIndex = spl.length > 1 ? spl.length - 1 : 1;
                   return resolveSecondary(spl.slice(undefined, endIndex).join("."), secPattern.slice(1));
                 }
                 return base + secPattern
               }
               return [
                       {
                           path: resolveSecondary(self.path, "^.bai"),
                           basename: resolveSecondary(self.basename, ".bai"),
                           class: "File",
                       }
               ];

       }
     outputBinding:
       glob: '$(inputs.outputBamName ? inputs.outputBamName : "generated")'
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - Mutect2
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 16, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4Mutect2


