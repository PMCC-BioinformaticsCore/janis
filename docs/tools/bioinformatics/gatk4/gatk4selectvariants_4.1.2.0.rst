:orphan:

GATK4: SelectVariants
===========================================

``Gatk4SelectVariants`` · *1 contributor · 4 versions*

USAGE: Selectvariants [arguments]
This tool makes it possible to select a subset of variants based on various criteria in order to facilitate certain
analyses. Examples include comparing and contrasting cases vs. controls, extracting variant or non-variant loci that
meet certain requirements, or troubleshooting some unexpected results, to name a few.
Version:4.1.3.0



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.selectvariants.versions import Gatk4SelectVariants_4_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4selectvariants_step",
           Gatk4SelectVariants_4_1_2(

           )
       )
       wf.output("out", source=gatk4selectvariants_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4SelectVariants:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4SelectVariants > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run Gatk4SelectVariants with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4SelectVariants





Information
------------

:ID: ``Gatk4SelectVariants``
:URL: *No URL to the documentation was provided*
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.2.0
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2019-12-02
:Updated: 2019-12-02


Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     Gzipped<VCF>
======  ============  ===============


Additional configuration (inputs)
---------------------------------

===================================  ==========================  =======================================  ==========  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                                 type                        prefix                                   position    documentation
===================================  ==========================  =======================================  ==========  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
javaOptions                          Optional<Array<String>>
compression_level                    Optional<Integer>                                                                Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
outputFilename                       Optional<Filename>          -O                                                   (--output) Required.
variants                             Optional<Gzipped<VCF>>      -V                                                   (--variant) A VCF file containing variants Required.
addOutputSamProgramRecord            Optional<Boolean>           -add-output-sam-program-record:Boolean               (--add-output-sam-program-record)  If true, adds a PG tag to created SAM/BAM/CRAM files.  Default value: true. Possible values: {true, false}
addOutputVcfCommandLine              Optional<Boolean>           -add-output-vcf-command-line                         (--add-output-vcf-command-line)  If true, adds a command line header line to created VCF files.  Default value: true. Possible values: {true, false}
arguments_file                       Optional<File>              --arguments_file                                     read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null.
cloudIndexPrefetchBuffer             Optional<Integer>           --cloud-index-prefetch-buffer                        (-CIPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.  Default value: -1.
cloudPrefetchBuffer                  Optional<Integer>           --cloud-prefetch-buffer                              (-CPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable).  Default value: 40.
conc                                 Optional<String>            -conc                                                (--concordance)  Output variants also called in this comparison track  Default value: null.
createOutputBamIndex                 Optional<Boolean>           --create-output-bam-index                            (-OBI)  If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.  Default value: true. Possible values: {true, false}
createOutputBamMd5                   Optional<Boolean>           --create-output-bam-md5                              (-OBM)  If true, create a MD5 digest for any BAM/SAM/CRAM file created  Default value: false. Possible values: {true, false}
createOutputvariantIndex             Optional<Boolean>           --create-output-variant-index                        (-OVI)  If true, create a VCF index when writing a coordinate-sorted VCF file.  Default value: true. Possible values: {true, false}
createOutputvariantMd5               Optional<Boolean>           --create-output-variant-md5                          (-OVM)  If true, create a a MD5 digest any VCF file created.  Default value: false. Possible values: {true, false}
disableBamIndexCaching               Optional<Boolean>           --disable-bam-index-caching                          (-DBIC:Boolean)  If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.  Default value: false. Possible values: {true, false}
disableReadFilter                    Optional<String>            --disable-read-filter                                (-DF)  Read filters to be disabled before analysis  This argument may be specified 0 or more times. Default value: null. Possible Values: {WellformedReadFilter}
disableSequenceDictionaryValidation  Optional<Boolean>           -disable-sequence-dictionary-validation              (--disable-sequence-dictionary-validation)  If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!  Default value: false. Possible values: {true, false}
disc                                 Optional<String>            -disc                                                (--discordance)  Output variants not called in this comparison track  Default value: null.
dropGenotypeAnnotation               Optional<String>            --drop-genotype-annotation                           (-DGA:String)  Genotype annotations to drop from output vcf.  Annotations to be dropped are specified by their key.  This argument may be specified 0 or more times. Default value: null.
dropInfoAnnotation                   Optional<String>            --drop-info-annotation                               (-DA:String)  Info annotations to drop from output vcf.  Annotations to be dropped are specified by their key.  This argument may be specified 0 or more times. Default value: null.
excludeFiltered                      Optional<Boolean>           --exclude-filtered                                   Don't include filtered sites Default value: false. Possible values: {true, false}
xlIds                                Optional<String>            -xl-ids                                              (--exclude-ids) List of variant rsIDs to exclude This argument may be specified 0 or more times. Default value: null.
excludeIntervals                     Optional<String>            --exclude-intervals                                  (-XL) This argument may be specified 0 or more times. Default value: null.
excludeNonvariants                   Optional<String>            --exclude-non-variants                               Default value: false. Possible values: {true, false}
excludeSampleExpressions             Optional<String>            --exclude-sample-expressions                         (-xl-se:String)  List of sample expressions to exclude  This argument may be specified 0 or more times. Default value: null.
excludeSampleName                    Optional<String>            --exclude-sample-name                                (-xl-sn:String)  Exclude genotypes from this sample  This argument may be specified 0 or more times. Default value: null.
gatkConfigFile                       Optional<File>              --gatk-config-file                                   A configuration file to use with the GATK. Default value: null.
gcsRetries                           Optional<Integer>           -gcs-retries                                         (--gcs-max-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20.
gcsProjectForRequesterPays           Optional<String>            --gcs-project-for-requester-pays                     Project to bill when accessing requester pays buckets. If unset, these buckets cannot be accessed.  Default value: .
help                                 Optional<Boolean>           -h                                                   (--help) display the help message Default value: false. Possible values: {true, false}
bam                                  Optional<IndexedBam>        -I                                                   (--input) BAM/SAM/CRAM file containing reads This argument may be specified 0 or more times. Default value: null.
intervalExclusionPadding             Optional<Integer>           --interval-exclusion-padding                         (-ixp:Integer)  Amount of padding (in bp) to add to each interval you are excluding.  Default value: 0.
imr                                  Optional<String>            -imr                                                 (--interval-merging-rule)  Interval merging rule for abutting intervals  Default value: ALL. Possible values: {ALL, OVERLAPPING_ONLY}
ip                                   Optional<Integer>           -ip                                                  (--interval-padding) Default value: 0.
isr                                  Optional<String>            -isr                                                 (--interval-set-rule)  Set merging approach to use for combining interval inputs  Default value: UNION. Possible values: {UNION, INTERSECTION}
intervals                            Optional<String>            --intervals                                          (-L:String) One or more genomic intervals over which to operate This argument may be specified 0 or more times. Default value: null.
invertMendelianViolation             Optional<Boolean>           --invert-mendelian-violation                         Output non-mendelian violation sites only  Default value: false. Possible values: {true, false}
invertSelect                         Optional<Boolean>           -invert-select                                       (--invertSelect)  Invert the selection criteria for -select  Default value: false. Possible values: {true, false}
ids                                  Optional<String>            -ids                                                 (--keep-ids) List of variant rsIDs to select This argument may be specified 0 or more times. Default value: null.
keepOriginalAc                       Optional<Boolean>           --keep-original-ac                                   Store the original AC, AF, and AN values after subsetting Default value: false. Possible values: {true, false}
keepOriginalDp                       Optional<Boolean>           --keep-original-dp                                   Store the original DP value after subsetting Default value: false. Possible values: {true, false}
le                                   Optional<Boolean>           -LE                                                  (--lenient) Lenient processing of VCF files Default value: false. Possible values: {true, false}
maxFilteredGenotypes                 Optional<Integer>           --max-filtered-genotypes                             Maximum number of samples filtered at the genotype level  Default value: 2147483647.
maxFractionFilteredGenotypes         Optional<Double>            --max-fraction-filtered-genotypes                    Maximum fraction of samples filtered at the genotype level  Default value: 1.0.
maxIndelSize                         Optional<Integer>           --max-indel-size                                     Maximum size of indels to include Default value: 2147483647.
maxNocallFraction                    Optional<Double>            --max-nocall-fraction                                Maximum fraction of samples with no-call genotypes Default value: 1.0.
maxNocallNumber                      Optional<Integer>           --max-nocall-number                                  Maximum number of samples with no-call genotypes Default value: 2147483647.
mendelianViolation                   Optional<Boolean>           --mendelian-violation                                Default value: false. Possible values: {true, false}
mendelianViolationQualThreshold      Optional<Double>            --mendelian-violation-qual-threshold                 Minimum GQ score for each trio member to accept a site as a violation  Default value: 0.0.
minFilteredGenotypes                 Optional<Integer>           --min-filtered-genotypes                             Minimum number of samples filtered at the genotype level  Default value: 0.
minFractionFilteredGenotypes         Optional<Double>            --min-fraction-filtered-genotypes                    Maximum fraction of samples filtered at the genotype level  Default value: 0.0.
minIndelSize                         Optional<Integer>           --min-indel-size                                     Minimum size of indels to include Default value: 0.
pedigree                             Optional<File>              --pedigree                                           (-ped:File) Pedigree file Default value: null.
preserveAlleles                      Optional<Boolean>           --preserve-alleles                                   Preserve original alleles, do not trim Default value: false. Possible values: {true, false}
quiet                                Optional<Boolean>           --QUIET                                              Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false}
readFilter                           Optional<String>            --read-filter                                        (-RF:String) Read filters to be applied before analysis This argument may be specified 0 or more times. Default value: null. Possible Values: {AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, SoftClippedReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter}
readIndex                            Optional<File>              -read-index                                          (--read-index)  Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically.  This argument may be specified 0 or more times. Default value: null.
readValidationStringency             Optional<String>            --read-validation-stringency                         (-VS:ValidationStringency)  Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: SILENT. Possible values: {STRICT, LENIENT, SILENT}
reference                            Optional<FastaWithIndexes>  --reference                                          (-R:String) Reference sequence Default value: null.
removeFractionGenotypes              Optional<Double>            --remove-fraction-genotypes                          Select a fraction of genotypes at random from the input and sets them to no-call  Default value: 0.0.
removeUnusedAlternates               Optional<Boolean>           --remove-unused-alternates                           Remove alternate alleles not present in any genotypes  Default value: false. Possible values: {true, false}
restrictAllelesTo                    Optional<String>            --restrict-alleles-to                                Select only variants of a particular allelicity  Default value: ALL. Possible values: {ALL, BIALLELIC, MULTIALLELIC}
sampleExpressions                    Optional<String>            --sample-expressions                                 (-se:String)  Regular expression to select multiple samples  This argument may be specified 0 or more times. Default value: null.
sampleName                           Optional<String>            --sample-name                                        (-sn:String) Include genotypes from this sample This argument may be specified 0 or more times. Default value: null.
secondsBetweenProgressUpdates        Optional<Double>            -seconds-between-progress-updates                    (--seconds-between-progress-updates)  Output traversal statistics every time this many seconds elapse  Default value: 10.0.
selectRandomFraction                 Optional<String>            --select-random-fraction                             (-fraction:Double)  Select a fraction of variants at random from the input  Default value: 0.0.
selectTypeToExclude                  Optional<String>            --select-type-to-exclude                             (-xl-select-type:Type)  Do not select certain type of variants from the input file  This argument may be specified 0 or more times. Default value: null. Possible values: {NO_VARIATION, SNP, MNP, INDEL, SYMBOLIC, MIXED}
selectTypeToInclude                  Optional<String>            --select-type-to-include                             (-select-type:Type)  Select only a certain type of variants from the input file  This argument may be specified 0 or more times. Default value: null. Possible values: {NO_VARIATION, SNP, MNP, INDEL, SYMBOLIC, MIXED}
selectexpressions                    Optional<String>            --selectExpressions                                  (-select:String)  One or more criteria to use when selecting the data  This argument may be specified 0 or more times. Default value: null.
sequenceDictionary                   Optional<File>              -sequence-dictionary                                 (--sequence-dictionary)  Use the given sequence dictionary as the master/canonical sequence dictionary.  Must be a .dict file.  Default value: null.
setFilteredGtToNocall                Optional<Boolean>           --set-filtered-gt-to-nocall                          Set filtered genotypes to no-call  Default value: false. Possible values: {true, false}
sitesOnlyVcfOutput                   Optional<Boolean>           --sites-only-vcf-output                              If true, don't emit genotype fields when writing vcf file output.  Default value: false. Possible values: {true, false}
tmpDir                               Optional<Filename>          --tmp-dir                                            Temp directory to use. Default value: null.
jdkDeflater                          Optional<Boolean>           -jdk-deflater                                        (--use-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false}
jdkInflater                          Optional<Boolean>           -jdk-inflater                                        (--use-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false}
verbosity                            Optional<String>            -verbosity                                           (--verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG}
version                              Optional<Boolean>           --version                                            display the version number for this tool Default value: false. Possible values: {true, false}
disableToolDefaultReadFilters        Optional<Boolean>           -disable-tool-default-read-filters                   (--disable-tool-default-read-filters)  Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)  Default value: false. Possible values: {true, false}
showhidden                           Optional<Boolean>           -showHidden                                          (--showHidden)  display hidden arguments  Default value: false. Possible values: {true, false}
ambigFilterBases                     Optional<Integer>           --ambig-filter-bases                                 Valid only if 'AmbiguousBaseReadFilter' is specified: Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction.  Default value: null.  Cannot be used in conjuction with argument(s) maxAmbiguousBaseFraction
ambigFilterFrac                      Optional<Double>            --ambig-filter-frac                                  Valid only if 'AmbiguousBaseReadFilter' is specified: Threshold fraction of ambiguous bases Default value: 0.05. Cannot be used in conjuction with argument(s) maxAmbiguousBases
maxFragmentLength                    Optional<Integer>           --max-fragment-length                                Valid only if 'FragmentLengthReadFilter' is specified: Maximum length of fragment (insert size) Default value: 1000000.
minFragmentLength                    Optional<Integer>           --min-fragment-length                                Valid only if 'FragmentLengthReadFilter' is specified: Minimum length of fragment (insert size) Default value: 0.
keepIntervals                        Optional<String>            --keep-intervals                                     Valid only if 'IntervalOverlapReadFilter' is specified: One or more genomic intervals to keep This argument must be specified at least once. Required.
library                              Optional<String>            -library                                             Valid only if 'LibraryReadFilter' is specified: (--library) Name of the library to keep This argument must be specified at least once. Required.
maximumMappingQuality                Optional<Integer>           --maximum-mapping-quality                            Valid only if 'MappingQualityReadFilter' is specified: Maximum mapping quality to keep (inclusive)  Default value: null.
minimumMappingQuality                Optional<Integer>           --minimum-mapping-quality                            Valid only if 'MappingQualityReadFilter' is specified: Minimum mapping quality to keep (inclusive)  Default value: 10.
dontRequireSoftClipsBothEnds         Optional<Boolean>           --dont-require-soft-clips-both-ends                  Valid only if 'OverclippedReadFilter' is specified: Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block  Default value: false. Possible values: {true, false}
filterTooShort                       Optional<Integer>           --filter-too-short                                   Valid only if 'OverclippedReadFilter' is specified: Minimum number of aligned bases Default value: 30.
platformFilterName                   Optional<String>            --platform-filter-name                               Valid only if 'PlatformReadFilter' is specified: This argument must be specified at least once. Required.
blackListedLanes                     Optional<String>            --black-listed-lanes                                 Valid only if 'PlatformUnitReadFilter' is specified: Platform unit (PU) to filter out This argument must be specified at least once. Required.
readGroupBlackList                   Optional<String>            --read-group-black-list                              Valid only if 'ReadGroupBlackListReadFilter' is specified: The name of the read group to filter out. This argument must be specified at least once. Required.
keepReadGroup                        Optional<String>            --keep-read-group                                    Valid only if 'ReadGroupReadFilter' is specified: The name of the read group to keep Required.
maxReadLength                        Optional<Integer>           --max-read-length                                    Valid only if 'ReadLengthReadFilter' is specified: Keep only reads with length at most equal to the specified value Required.
minReadLength                        Optional<Integer>           --min-read-length                                    Valid only if 'ReadLengthReadFilter' is specified: Keep only reads with length at least equal to the specified value Default value: 1.
readName                             Optional<String>            --read-name                                          Valid only if 'ReadNameReadFilter' is specified: Keep only reads with this read name Required.
keepReverseStrandOnly                Optional<Boolean>           --keep-reverse-strand-only                           Valid only if 'ReadStrandFilter' is specified: Keep only reads on the reverse strand  Required. Possible values: {true, false}
sample                               Optional<String>            --sample                                             Valid only if 'SampleReadFilter' is specified: The name of the sample(s) to keep, filtering out all others This argument must be specified at least once. Required.
invertSoftClipRatioFilter            Optional<Boolean>           --invert-soft-clip-ratio-filter                      Inverts the results from this filter, causing all variants that would pass to fail and visa-versa.  Default value: false. Possible values: {true, false}
softClippedLeadingTrailingRatio      Optional<Double>            --soft-clipped-leading-trailing-ratio                Threshold ratio of soft clipped bases (leading / trailing the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumSoftClippedRatio
softClippedRatioThreshold            Optional<Double>            --soft-clipped-ratio-threshold                       Threshold ratio of soft clipped bases (anywhere in the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumLeadingTrailingSoftClippedRatio
===================================  ==========================  =======================================  ==========  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4SelectVariants {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       String? outputFilename
       File? variants
       File? variants_tbi
       Boolean? addOutputSamProgramRecord
       Boolean? addOutputVcfCommandLine
       File? arguments_file
       Int? cloudIndexPrefetchBuffer
       Int? cloudPrefetchBuffer
       String? conc
       Boolean? createOutputBamIndex
       Boolean? createOutputBamMd5
       Boolean? createOutputvariantIndex
       Boolean? createOutputvariantMd5
       Boolean? disableBamIndexCaching
       String? disableReadFilter
       Boolean? disableSequenceDictionaryValidation
       String? disc
       String? dropGenotypeAnnotation
       String? dropInfoAnnotation
       Boolean? excludeFiltered
       String? xlIds
       String? excludeIntervals
       String? excludeNonvariants
       String? excludeSampleExpressions
       String? excludeSampleName
       File? gatkConfigFile
       Int? gcsRetries
       String? gcsProjectForRequesterPays
       Boolean? help
       File? bam
       File? bam_bai
       Int? intervalExclusionPadding
       String? imr
       Int? ip
       String? isr
       String? intervals
       Boolean? invertMendelianViolation
       Boolean? invertSelect
       String? ids
       Boolean? keepOriginalAc
       Boolean? keepOriginalDp
       Boolean? le
       Int? maxFilteredGenotypes
       Float? maxFractionFilteredGenotypes
       Int? maxIndelSize
       Float? maxNocallFraction
       Int? maxNocallNumber
       Boolean? mendelianViolation
       Float? mendelianViolationQualThreshold
       Int? minFilteredGenotypes
       Float? minFractionFilteredGenotypes
       Int? minIndelSize
       File? pedigree
       Boolean? preserveAlleles
       Boolean? quiet
       String? readFilter
       File? readIndex
       String? readValidationStringency
       File? reference
       File? reference_fai
       File? reference_amb
       File? reference_ann
       File? reference_bwt
       File? reference_pac
       File? reference_sa
       File? reference_dict
       Float? removeFractionGenotypes
       Boolean? removeUnusedAlternates
       String? restrictAllelesTo
       String? sampleExpressions
       String? sampleName
       Float? secondsBetweenProgressUpdates
       String? selectRandomFraction
       String? selectTypeToExclude
       String? selectTypeToInclude
       String? selectexpressions
       File? sequenceDictionary
       Boolean? setFilteredGtToNocall
       Boolean? sitesOnlyVcfOutput
       String? tmpDir
       Boolean? jdkDeflater
       Boolean? jdkInflater
       String? verbosity
       Boolean? version
       Boolean? disableToolDefaultReadFilters
       Boolean? showhidden
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
       Boolean? invertSoftClipRatioFilter
       Float? softClippedLeadingTrailingRatio
       Float? softClippedRatioThreshold
     }
     command <<<
       set -e
       gatk SelectVariants \
         --java-options '-Xmx~{((select_first([runtime_memory, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         -O '~{select_first([outputFilename, "generated"])}' \
         ~{if defined(variants) then ("-V '" + variants + "'") else ""} \
         ~{if (defined(addOutputSamProgramRecord) && select_first([addOutputSamProgramRecord])) then "-add-output-sam-program-record:Boolean" else ""} \
         ~{if (defined(addOutputVcfCommandLine) && select_first([addOutputVcfCommandLine])) then "-add-output-vcf-command-line" else ""} \
         ~{if defined(arguments_file) then ("--arguments_file '" + arguments_file + "'") else ""} \
         ~{if defined(cloudIndexPrefetchBuffer) then ("--cloud-index-prefetch-buffer " + cloudIndexPrefetchBuffer) else ''} \
         ~{if defined(cloudPrefetchBuffer) then ("--cloud-prefetch-buffer " + cloudPrefetchBuffer) else ''} \
         ~{if defined(conc) then ("-conc '" + conc + "'") else ""} \
         ~{if (defined(createOutputBamIndex) && select_first([createOutputBamIndex])) then "--create-output-bam-index" else ""} \
         ~{if (defined(createOutputBamMd5) && select_first([createOutputBamMd5])) then "--create-output-bam-md5" else ""} \
         ~{if select_first([createOutputvariantIndex, true]) then "--create-output-variant-index" else ""} \
         ~{if (defined(createOutputvariantMd5) && select_first([createOutputvariantMd5])) then "--create-output-variant-md5" else ""} \
         ~{if (defined(disableBamIndexCaching) && select_first([disableBamIndexCaching])) then "--disable-bam-index-caching" else ""} \
         ~{if defined(disableReadFilter) then ("--disable-read-filter '" + disableReadFilter + "'") else ""} \
         ~{if (defined(disableSequenceDictionaryValidation) && select_first([disableSequenceDictionaryValidation])) then "-disable-sequence-dictionary-validation" else ""} \
         ~{if defined(disc) then ("-disc '" + disc + "'") else ""} \
         ~{if defined(dropGenotypeAnnotation) then ("--drop-genotype-annotation '" + dropGenotypeAnnotation + "'") else ""} \
         ~{if defined(dropInfoAnnotation) then ("--drop-info-annotation '" + dropInfoAnnotation + "'") else ""} \
         ~{if (defined(excludeFiltered) && select_first([excludeFiltered])) then "--exclude-filtered" else ""} \
         ~{if defined(xlIds) then ("-xl-ids '" + xlIds + "'") else ""} \
         ~{if defined(excludeIntervals) then ("--exclude-intervals '" + excludeIntervals + "'") else ""} \
         ~{if defined(excludeNonvariants) then ("--exclude-non-variants '" + excludeNonvariants + "'") else ""} \
         ~{if defined(excludeSampleExpressions) then ("--exclude-sample-expressions '" + excludeSampleExpressions + "'") else ""} \
         ~{if defined(excludeSampleName) then ("--exclude-sample-name '" + excludeSampleName + "'") else ""} \
         ~{if defined(gatkConfigFile) then ("--gatk-config-file '" + gatkConfigFile + "'") else ""} \
         ~{if defined(gcsRetries) then ("-gcs-retries " + gcsRetries) else ''} \
         ~{if defined(gcsProjectForRequesterPays) then ("--gcs-project-for-requester-pays '" + gcsProjectForRequesterPays + "'") else ""} \
         ~{if (defined(help) && select_first([help])) then "-h" else ""} \
         ~{if defined(bam) then ("-I '" + bam + "'") else ""} \
         ~{if defined(intervalExclusionPadding) then ("--interval-exclusion-padding " + intervalExclusionPadding) else ''} \
         ~{if defined(imr) then ("-imr '" + imr + "'") else ""} \
         ~{if defined(ip) then ("-ip " + ip) else ''} \
         ~{if defined(isr) then ("-isr '" + isr + "'") else ""} \
         ~{if defined(intervals) then ("--intervals '" + intervals + "'") else ""} \
         ~{if (defined(invertMendelianViolation) && select_first([invertMendelianViolation])) then "--invert-mendelian-violation" else ""} \
         ~{if (defined(invertSelect) && select_first([invertSelect])) then "-invert-select" else ""} \
         ~{if defined(ids) then ("-ids '" + ids + "'") else ""} \
         ~{if (defined(keepOriginalAc) && select_first([keepOriginalAc])) then "--keep-original-ac" else ""} \
         ~{if (defined(keepOriginalDp) && select_first([keepOriginalDp])) then "--keep-original-dp" else ""} \
         ~{if (defined(le) && select_first([le])) then "-LE" else ""} \
         ~{if defined(maxFilteredGenotypes) then ("--max-filtered-genotypes " + maxFilteredGenotypes) else ''} \
         ~{if defined(maxFractionFilteredGenotypes) then ("--max-fraction-filtered-genotypes " + maxFractionFilteredGenotypes) else ''} \
         ~{if defined(maxIndelSize) then ("--max-indel-size " + maxIndelSize) else ''} \
         ~{if defined(maxNocallFraction) then ("--max-nocall-fraction " + maxNocallFraction) else ''} \
         ~{if defined(maxNocallNumber) then ("--max-nocall-number " + maxNocallNumber) else ''} \
         ~{if (defined(mendelianViolation) && select_first([mendelianViolation])) then "--mendelian-violation" else ""} \
         ~{if defined(mendelianViolationQualThreshold) then ("--mendelian-violation-qual-threshold " + mendelianViolationQualThreshold) else ''} \
         ~{if defined(minFilteredGenotypes) then ("--min-filtered-genotypes " + minFilteredGenotypes) else ''} \
         ~{if defined(minFractionFilteredGenotypes) then ("--min-fraction-filtered-genotypes " + minFractionFilteredGenotypes) else ''} \
         ~{if defined(minIndelSize) then ("--min-indel-size " + minIndelSize) else ''} \
         ~{if defined(pedigree) then ("--pedigree '" + pedigree + "'") else ""} \
         ~{if (defined(preserveAlleles) && select_first([preserveAlleles])) then "--preserve-alleles" else ""} \
         ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
         ~{if defined(readFilter) then ("--read-filter '" + readFilter + "'") else ""} \
         ~{if defined(readIndex) then ("-read-index '" + readIndex + "'") else ""} \
         ~{if defined(readValidationStringency) then ("--read-validation-stringency '" + readValidationStringency + "'") else ""} \
         ~{if defined(reference) then ("--reference '" + reference + "'") else ""} \
         ~{if defined(removeFractionGenotypes) then ("--remove-fraction-genotypes " + removeFractionGenotypes) else ''} \
         ~{if (defined(removeUnusedAlternates) && select_first([removeUnusedAlternates])) then "--remove-unused-alternates" else ""} \
         ~{if defined(restrictAllelesTo) then ("--restrict-alleles-to '" + restrictAllelesTo + "'") else ""} \
         ~{if defined(sampleExpressions) then ("--sample-expressions '" + sampleExpressions + "'") else ""} \
         ~{if defined(sampleName) then ("--sample-name '" + sampleName + "'") else ""} \
         ~{if defined(secondsBetweenProgressUpdates) then ("-seconds-between-progress-updates " + secondsBetweenProgressUpdates) else ''} \
         ~{if defined(selectRandomFraction) then ("--select-random-fraction '" + selectRandomFraction + "'") else ""} \
         ~{if defined(selectTypeToExclude) then ("--select-type-to-exclude '" + selectTypeToExclude + "'") else ""} \
         ~{if defined(selectTypeToInclude) then ("--select-type-to-include '" + selectTypeToInclude + "'") else ""} \
         ~{if defined(selectexpressions) then ("--selectExpressions '" + selectexpressions + "'") else ""} \
         ~{if defined(sequenceDictionary) then ("-sequence-dictionary '" + sequenceDictionary + "'") else ""} \
         ~{if (defined(setFilteredGtToNocall) && select_first([setFilteredGtToNocall])) then "--set-filtered-gt-to-nocall" else ""} \
         ~{if (defined(sitesOnlyVcfOutput) && select_first([sitesOnlyVcfOutput])) then "--sites-only-vcf-output" else ""} \
         --tmp-dir '~{select_first([tmpDir, "generated"])}' \
         ~{if (defined(jdkDeflater) && select_first([jdkDeflater])) then "-jdk-deflater" else ""} \
         ~{if (defined(jdkInflater) && select_first([jdkInflater])) then "-jdk-inflater" else ""} \
         ~{if defined(verbosity) then ("-verbosity '" + verbosity + "'") else ""} \
         ~{if (defined(version) && select_first([version])) then "--version" else ""} \
         ~{if (defined(disableToolDefaultReadFilters) && select_first([disableToolDefaultReadFilters])) then "-disable-tool-default-read-filters" else ""} \
         ~{if (defined(showhidden) && select_first([showhidden])) then "-showHidden" else ""} \
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
         ~{if defined(readGroupBlackList) then ("--read-group-black-list '" + readGroupBlackList + "'") else ""} \
         ~{if defined(keepReadGroup) then ("--keep-read-group '" + keepReadGroup + "'") else ""} \
         ~{if defined(maxReadLength) then ("--max-read-length " + maxReadLength) else ''} \
         ~{if defined(minReadLength) then ("--min-read-length " + minReadLength) else ''} \
         ~{if defined(readName) then ("--read-name '" + readName + "'") else ""} \
         ~{if (defined(keepReverseStrandOnly) && select_first([keepReverseStrandOnly])) then "--keep-reverse-strand-only" else ""} \
         ~{if defined(sample) then ("--sample '" + sample + "'") else ""} \
         ~{if (defined(invertSoftClipRatioFilter) && select_first([invertSoftClipRatioFilter])) then "--invert-soft-clip-ratio-filter" else ""} \
         ~{if defined(softClippedLeadingTrailingRatio) then ("--soft-clipped-leading-trailing-ratio " + softClippedLeadingTrailingRatio) else ''} \
         ~{if defined(softClippedRatioThreshold) then ("--soft-clipped-ratio-threshold " + softClippedRatioThreshold) else ''}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.2.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated"])
       File out_tbi = select_first([outputFilename, "generated"]) + ".tbi"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: SelectVariants'
   doc: |
     USAGE: Selectvariants [arguments]
     This tool makes it possible to select a subset of variants based on various criteria in order to facilitate certain
     analyses. Examples include comparing and contrasting cases vs. controls, extracting variant or non-variant loci that
     meet certain requirements, or troubleshooting some unexpected results, to name a few.
     Version:4.1.3.0

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.2.0

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
   - id: outputFilename
     label: outputFilename
     doc: (--output) Required.
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: -O
       separate: true
   - id: variants
     label: variants
     doc: (--variant) A VCF file containing variants Required.
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       prefix: -V
       separate: true
   - id: addOutputSamProgramRecord
     label: addOutputSamProgramRecord
     doc: |-
       (--add-output-sam-program-record)  If true, adds a PG tag to created SAM/BAM/CRAM files.  Default value: true. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -add-output-sam-program-record:Boolean
       separate: true
   - id: addOutputVcfCommandLine
     label: addOutputVcfCommandLine
     doc: |-
       (--add-output-vcf-command-line)  If true, adds a command line header line to created VCF files.  Default value: true. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -add-output-vcf-command-line
       separate: true
   - id: arguments_file
     label: arguments_file
     doc: |-
       read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null. 
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --arguments_file
       separate: true
   - id: cloudIndexPrefetchBuffer
     label: cloudIndexPrefetchBuffer
     doc: |-
       (-CIPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.  Default value: -1. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --cloud-index-prefetch-buffer
       separate: true
   - id: cloudPrefetchBuffer
     label: cloudPrefetchBuffer
     doc: |-
       (-CPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable).  Default value: 40. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --cloud-prefetch-buffer
       separate: true
   - id: conc
     label: conc
     doc: |-
       (--concordance)  Output variants also called in this comparison track  Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -conc
       separate: true
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
   - id: createOutputvariantIndex
     label: createOutputvariantIndex
     doc: |-
       (-OVI)  If true, create a VCF index when writing a coordinate-sorted VCF file.  Default value: true. Possible values: {true, false} 
     type: boolean
     default: true
     inputBinding:
       prefix: --create-output-variant-index
   - id: createOutputvariantMd5
     label: createOutputvariantMd5
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
       (-DBIC:Boolean)  If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --disable-bam-index-caching
   - id: disableReadFilter
     label: disableReadFilter
     doc: |-
       (-DF)  Read filters to be disabled before analysis  This argument may be specified 0 or more times. Default value: null. Possible Values: {WellformedReadFilter}
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --disable-read-filter
       separate: true
   - id: disableSequenceDictionaryValidation
     label: disableSequenceDictionaryValidation
     doc: |-
       (--disable-sequence-dictionary-validation)  If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -disable-sequence-dictionary-validation
   - id: disc
     label: disc
     doc: |-
       (--discordance)  Output variants not called in this comparison track  Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -disc
       separate: true
   - id: dropGenotypeAnnotation
     label: dropGenotypeAnnotation
     doc: |-
       (-DGA:String)  Genotype annotations to drop from output vcf.  Annotations to be dropped are specified by their key.  This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --drop-genotype-annotation
       separate: true
   - id: dropInfoAnnotation
     label: dropInfoAnnotation
     doc: |-
       (-DA:String)  Info annotations to drop from output vcf.  Annotations to be dropped are specified by their key.  This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --drop-info-annotation
       separate: true
   - id: excludeFiltered
     label: excludeFiltered
     doc: |-
       Don't include filtered sites Default value: false. Possible values: {true, false}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --exclude-filtered
   - id: xlIds
     label: xlIds
     doc: |-
       (--exclude-ids) List of variant rsIDs to exclude This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -xl-ids
       separate: true
   - id: excludeIntervals
     label: excludeIntervals
     doc: '(-XL) This argument may be specified 0 or more times. Default value: null. '
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --exclude-intervals
       separate: true
   - id: excludeNonvariants
     label: excludeNonvariants
     doc: 'Default value: false. Possible values: {true, false}'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --exclude-non-variants
   - id: excludeSampleExpressions
     label: excludeSampleExpressions
     doc: |-
       (-xl-se:String)  List of sample expressions to exclude  This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --exclude-sample-expressions
       separate: true
   - id: excludeSampleName
     label: excludeSampleName
     doc: |-
       (-xl-sn:String)  Exclude genotypes from this sample  This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --exclude-sample-name
       separate: true
   - id: gatkConfigFile
     label: gatkConfigFile
     doc: 'A configuration file to use with the GATK. Default value: null.'
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --gatk-config-file
       separate: true
   - id: gcsRetries
     label: gcsRetries
     doc: |-
       (--gcs-max-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -gcs-retries
       separate: true
   - id: gcsProjectForRequesterPays
     label: gcsProjectForRequesterPays
     doc: |2-
        Project to bill when accessing requester pays buckets. If unset, these buckets cannot be accessed.  Default value: . 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --gcs-project-for-requester-pays
       separate: true
   - id: help
     label: help
     doc: |-
       (--help) display the help message Default value: false. Possible values: {true, false}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -h
   - id: bam
     label: bam
     doc: |-
       (--input) BAM/SAM/CRAM file containing reads This argument may be specified 0 or more times. Default value: null. 
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       prefix: -I
       separate: true
   - id: intervalExclusionPadding
     label: intervalExclusionPadding
     doc: |-
       (-ixp:Integer)  Amount of padding (in bp) to add to each interval you are excluding.  Default value: 0. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --interval-exclusion-padding
       separate: true
   - id: imr
     label: imr
     doc: |-
       (--interval-merging-rule)  Interval merging rule for abutting intervals  Default value: ALL. Possible values: {ALL, OVERLAPPING_ONLY} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -imr
       separate: true
   - id: ip
     label: ip
     doc: '(--interval-padding) Default value: 0.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -ip
       separate: true
   - id: isr
     label: isr
     doc: |-
       (--interval-set-rule)  Set merging approach to use for combining interval inputs  Default value: UNION. Possible values: {UNION, INTERSECTION} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -isr
       separate: true
   - id: intervals
     label: intervals
     doc: |-
       (-L:String) One or more genomic intervals over which to operate This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --intervals
       separate: true
   - id: invertMendelianViolation
     label: invertMendelianViolation
     doc: |2-
        Output non-mendelian violation sites only  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --invert-mendelian-violation
   - id: invertSelect
     label: invertSelect
     doc: |-
       (--invertSelect)  Invert the selection criteria for -select  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -invert-select
   - id: ids
     label: ids
     doc: |-
       (--keep-ids) List of variant rsIDs to select This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -ids
       separate: true
   - id: keepOriginalAc
     label: keepOriginalAc
     doc: |-
       Store the original AC, AF, and AN values after subsetting Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --keep-original-ac
   - id: keepOriginalDp
     label: keepOriginalDp
     doc: |-
       Store the original DP value after subsetting Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --keep-original-dp
       separate: true
   - id: le
     label: le
     doc: |-
       (--lenient) Lenient processing of VCF files Default value: false. Possible values: {true, false}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -LE
       separate: true
   - id: maxFilteredGenotypes
     label: maxFilteredGenotypes
     doc: |-
       Maximum number of samples filtered at the genotype level  Default value: 2147483647. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-filtered-genotypes
       separate: true
   - id: maxFractionFilteredGenotypes
     label: maxFractionFilteredGenotypes
     doc: |2-
        Maximum fraction of samples filtered at the genotype level  Default value: 1.0. 
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --max-fraction-filtered-genotypes
       separate: true
   - id: maxIndelSize
     label: maxIndelSize
     doc: 'Maximum size of indels to include Default value: 2147483647.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-indel-size
       separate: true
   - id: maxNocallFraction
     label: maxNocallFraction
     doc: 'Maximum fraction of samples with no-call genotypes Default value: 1.0.'
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --max-nocall-fraction
       separate: true
   - id: maxNocallNumber
     label: maxNocallNumber
     doc: 'Maximum number of samples with no-call genotypes Default value: 2147483647.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-nocall-number
       separate: true
   - id: mendelianViolation
     label: mendelianViolation
     doc: 'Default value: false. Possible values: {true, false} '
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --mendelian-violation
       separate: true
   - id: mendelianViolationQualThreshold
     label: mendelianViolationQualThreshold
     doc: |2-
        Minimum GQ score for each trio member to accept a site as a violation  Default value: 0.0.
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --mendelian-violation-qual-threshold
       separate: true
   - id: minFilteredGenotypes
     label: minFilteredGenotypes
     doc: ' Minimum number of samples filtered at the genotype level  Default value:
       0. '
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-filtered-genotypes
       separate: true
   - id: minFractionFilteredGenotypes
     label: minFractionFilteredGenotypes
     doc: |2-
        Maximum fraction of samples filtered at the genotype level  Default value: 0.0. 
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --min-fraction-filtered-genotypes
       separate: true
   - id: minIndelSize
     label: minIndelSize
     doc: 'Minimum size of indels to include Default value: 0.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-indel-size
       separate: true
   - id: pedigree
     label: pedigree
     doc: '(-ped:File) Pedigree file Default value: null.'
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --pedigree
       separate: true
   - id: preserveAlleles
     label: preserveAlleles
     doc: |-
       Preserve original alleles, do not trim Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --preserve-alleles
       separate: true
   - id: quiet
     label: quiet
     doc: |-
       Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --QUIET
       separate: true
   - id: readFilter
     label: readFilter
     doc: |-
       (-RF:String) Read filters to be applied before analysis This argument may be specified 0 or more times. Default value: null. Possible Values: {AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, SoftClippedReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter}
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read-filter
       separate: true
   - id: readIndex
     label: readIndex
     doc: |-
       (--read-index)  Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically.  This argument may be specified 0 or more times. Default value: null. 
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -read-index
       separate: true
   - id: readValidationStringency
     label: readValidationStringency
     doc: |-
       (-VS:ValidationStringency)  Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: SILENT. Possible values: {STRICT, LENIENT, SILENT} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read-validation-stringency
       separate: true
   - id: reference
     label: reference
     doc: '(-R:String) Reference sequence Default value: null.'
     type:
     - File
     - 'null'
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
       separate: true
   - id: removeFractionGenotypes
     label: removeFractionGenotypes
     doc: |2-
        Select a fraction of genotypes at random from the input and sets them to no-call  Default value: 0.0. 
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --remove-fraction-genotypes
       separate: true
   - id: removeUnusedAlternates
     label: removeUnusedAlternates
     doc: |2-
        Remove alternate alleles not present in any genotypes  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --remove-unused-alternates
       separate: true
   - id: restrictAllelesTo
     label: restrictAllelesTo
     doc: |2-
        Select only variants of a particular allelicity  Default value: ALL. Possible values: {ALL, BIALLELIC, MULTIALLELIC} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --restrict-alleles-to
       separate: true
   - id: sampleExpressions
     label: sampleExpressions
     doc: |-
       (-se:String)  Regular expression to select multiple samples  This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sample-expressions
       separate: true
   - id: sampleName
     label: sampleName
     doc: |-
       (-sn:String) Include genotypes from this sample This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sample-name
       separate: true
   - id: secondsBetweenProgressUpdates
     label: secondsBetweenProgressUpdates
     doc: |-
       (--seconds-between-progress-updates)  Output traversal statistics every time this many seconds elapse  Default value: 10.0. 
     type:
     - double
     - 'null'
     inputBinding:
       prefix: -seconds-between-progress-updates
       separate: true
   - id: selectRandomFraction
     label: selectRandomFraction
     doc: |-
       (-fraction:Double)  Select a fraction of variants at random from the input  Default value: 0.0. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --select-random-fraction
       separate: true
   - id: selectTypeToExclude
     label: selectTypeToExclude
     doc: |-
       (-xl-select-type:Type)  Do not select certain type of variants from the input file  This argument may be specified 0 or more times. Default value: null. Possible values: {NO_VARIATION, SNP, MNP, INDEL, SYMBOLIC, MIXED} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --select-type-to-exclude
       separate: true
   - id: selectTypeToInclude
     label: selectTypeToInclude
     doc: |-
       (-select-type:Type)  Select only a certain type of variants from the input file  This argument may be specified 0 or more times. Default value: null. Possible values: {NO_VARIATION, SNP, MNP, INDEL, SYMBOLIC, MIXED} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --select-type-to-include
       separate: true
   - id: selectexpressions
     label: selectexpressions
     doc: |-
       (-select:String)  One or more criteria to use when selecting the data  This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --selectExpressions
       separate: true
   - id: sequenceDictionary
     label: sequenceDictionary
     doc: |-
       (--sequence-dictionary)  Use the given sequence dictionary as the master/canonical sequence dictionary.  Must be a .dict file.  Default value: null. 
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -sequence-dictionary
       separate: true
   - id: setFilteredGtToNocall
     label: setFilteredGtToNocall
     doc: |2-
        Set filtered genotypes to no-call  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --set-filtered-gt-to-nocall
       separate: true
   - id: sitesOnlyVcfOutput
     label: sitesOnlyVcfOutput
     doc: |2-
        If true, don't emit genotype fields when writing vcf file output.  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --sites-only-vcf-output
       separate: true
   - id: tmpDir
     label: tmpDir
     doc: 'Temp directory to use. Default value: null.'
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --tmp-dir
       separate: true
   - id: jdkDeflater
     label: jdkDeflater
     doc: |-
       (--use-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -jdk-deflater
       separate: true
   - id: jdkInflater
     label: jdkInflater
     doc: |-
       (--use-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -jdk-inflater
       separate: true
   - id: verbosity
     label: verbosity
     doc: |-
       (--verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -verbosity
       separate: true
   - id: version
     label: version
     doc: |-
       display the version number for this tool Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --version
       separate: true
   - id: disableToolDefaultReadFilters
     label: disableToolDefaultReadFilters
     doc: |-
       (--disable-tool-default-read-filters)  Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -disable-tool-default-read-filters
       separate: true
   - id: showhidden
     label: showhidden
     doc: |-
       (--showHidden)  display hidden arguments  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -showHidden
       separate: true
   - id: ambigFilterBases
     label: ambigFilterBases
     doc: |-
       Valid only if 'AmbiguousBaseReadFilter' is specified: Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction.  Default value: null.  Cannot be used in conjuction with argument(s) maxAmbiguousBaseFraction
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --ambig-filter-bases
       separate: true
   - id: ambigFilterFrac
     label: ambigFilterFrac
     doc: |-
       Valid only if 'AmbiguousBaseReadFilter' is specified: Threshold fraction of ambiguous bases Default value: 0.05. Cannot be used in conjuction with argument(s) maxAmbiguousBases
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --ambig-filter-frac
       separate: true
   - id: maxFragmentLength
     label: maxFragmentLength
     doc: |-
       Valid only if 'FragmentLengthReadFilter' is specified: Maximum length of fragment (insert size) Default value: 1000000.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-fragment-length
       separate: true
   - id: minFragmentLength
     label: minFragmentLength
     doc: |-
       Valid only if 'FragmentLengthReadFilter' is specified: Minimum length of fragment (insert size) Default value: 0.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-fragment-length
       separate: true
   - id: keepIntervals
     label: keepIntervals
     doc: |-
       Valid only if 'IntervalOverlapReadFilter' is specified: One or more genomic intervals to keep This argument must be specified at least once. Required. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --keep-intervals
       separate: true
   - id: library
     label: library
     doc: |-
       Valid only if 'LibraryReadFilter' is specified: (--library) Name of the library to keep This argument must be specified at least once. Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -library
       separate: true
   - id: maximumMappingQuality
     label: maximumMappingQuality
     doc: |-
       Valid only if 'MappingQualityReadFilter' is specified: Maximum mapping quality to keep (inclusive)  Default value: null. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --maximum-mapping-quality
       separate: true
   - id: minimumMappingQuality
     label: minimumMappingQuality
     doc: |-
       Valid only if 'MappingQualityReadFilter' is specified: Minimum mapping quality to keep (inclusive)  Default value: 10. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --minimum-mapping-quality
       separate: true
   - id: dontRequireSoftClipsBothEnds
     label: dontRequireSoftClipsBothEnds
     doc: |-
       Valid only if 'OverclippedReadFilter' is specified: Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --dont-require-soft-clips-both-ends
       separate: true
   - id: filterTooShort
     label: filterTooShort
     doc: |-
       Valid only if 'OverclippedReadFilter' is specified: Minimum number of aligned bases Default value: 30.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --filter-too-short
       separate: true
   - id: platformFilterName
     label: platformFilterName
     doc: |-
       Valid only if 'PlatformReadFilter' is specified: This argument must be specified at least once. Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --platform-filter-name
       separate: true
   - id: blackListedLanes
     label: blackListedLanes
     doc: |-
       Valid only if 'PlatformUnitReadFilter' is specified: Platform unit (PU) to filter out This argument must be specified at least once. Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --black-listed-lanes
       separate: true
   - id: readGroupBlackList
     label: readGroupBlackList
     doc: |-
       Valid only if 'ReadGroupBlackListReadFilter' is specified: The name of the read group to filter out. This argument must be specified at least once. Required. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read-group-black-list
       separate: true
   - id: keepReadGroup
     label: keepReadGroup
     doc: |-
       Valid only if 'ReadGroupReadFilter' is specified: The name of the read group to keep Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --keep-read-group
       separate: true
   - id: maxReadLength
     label: maxReadLength
     doc: |-
       Valid only if 'ReadLengthReadFilter' is specified: Keep only reads with length at most equal to the specified value Required.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-read-length
       separate: true
   - id: minReadLength
     label: minReadLength
     doc: |-
       Valid only if 'ReadLengthReadFilter' is specified: Keep only reads with length at least equal to the specified value Default value: 1.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-read-length
       separate: true
   - id: readName
     label: readName
     doc: |-
       Valid only if 'ReadNameReadFilter' is specified: Keep only reads with this read name Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read-name
       separate: true
   - id: keepReverseStrandOnly
     label: keepReverseStrandOnly
     doc: |-
       Valid only if 'ReadStrandFilter' is specified: Keep only reads on the reverse strand  Required. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --keep-reverse-strand-only
       separate: true
   - id: sample
     label: sample
     doc: |-
       Valid only if 'SampleReadFilter' is specified: The name of the sample(s) to keep, filtering out all others This argument must be specified at least once. Required. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sample
       separate: true
   - id: invertSoftClipRatioFilter
     label: invertSoftClipRatioFilter
     doc: |2-
        Inverts the results from this filter, causing all variants that would pass to fail and visa-versa.  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --invert-soft-clip-ratio-filter
       separate: true
   - id: softClippedLeadingTrailingRatio
     label: softClippedLeadingTrailingRatio
     doc: |2-
        Threshold ratio of soft clipped bases (leading / trailing the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumSoftClippedRatio
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --soft-clipped-leading-trailing-ratio
       separate: true
   - id: softClippedRatioThreshold
     label: softClippedRatioThreshold
     doc: |2-
        Threshold ratio of soft clipped bases (anywhere in the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumLeadingTrailingSoftClippedRatio
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --soft-clipped-ratio-threshold
       separate: true

   outputs:
   - id: out
     label: out
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: generated
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - SelectVariants
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4SelectVariants


