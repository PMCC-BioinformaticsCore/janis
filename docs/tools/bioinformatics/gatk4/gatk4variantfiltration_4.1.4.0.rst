:orphan:

GATK4: VariantFiltration
=================================================

``Gatk4Variantfiltration`` · *1 contributor · 4 versions*

USAGE: VariantFiltration [arguments]
Filter variant calls based on INFO and/or FORMAT annotations.
Version:4.1.3.0



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.variantfiltration.versions import Gatk4VariantFiltration_4_1_4

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4variantfiltration_step",
           Gatk4VariantFiltration_4_1_4(
               variant=None,
               filterName=None,
           )
       )
       wf.output("out", source=gatk4variantfiltration_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4Variantfiltration:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4Variantfiltration > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       filterName:
       - <value>
       - <value>
       variant: variant.vcf.gz




5. Run Gatk4Variantfiltration with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4Variantfiltration





Information
------------

:ID: ``Gatk4Variantfiltration``
:URL: *No URL to the documentation was provided*
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.4.0
:Authors: Michael Franklin
:Citations: None
:Created: 2020-05-18
:Updated: 2020-05-18


Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     Gzipped<VCF>
======  ============  ===============


Additional configuration (inputs)
---------------------------------

===================================  ==========================  ========================================  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                                 type                        prefix                                    position    documentation
===================================  ==========================  ========================================  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
variant                              Gzipped<VCF>                --variant                                             (-V) A VCF file containing variants Required.
filterName                           Array<Optional<String>>     --filter-name                                         Names to use for the list of filters This argument may be specified 0 or more times. Default value: null.
javaOptions                          Optional<Array<String>>
compression_level                    Optional<Integer>                                                                 Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
outputFilename                       Optional<Filename>          --output                                              (-O) File to which variants should be written Required.
addOutputSamProgramRecord            Optional<Boolean>           --add-output-sam-program-record                       (-add-output-sam-program-record)  If true, adds a PG tag to created SAM/BAM/CRAM files.  Default value: true. Possible values: {true, false}
addOutputVcfCommandLine              Optional<Boolean>           --add-output-vcf-command-line                         (-add-output-vcf-command-line)  If true, adds a command line header line to created VCF files.  Default value: true. Possible values: {true, false}
arguments_file                       Optional<File>              --arguments_file                                      read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null.
cloudIndexPrefetchBuffer             Optional<Integer>           --cloud-index-prefetch-buffer                         (-CIPB)  Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.  Default value: -1.
cloudPrefetchBuffer                  Optional<Integer>           --cloud-prefetch-buffer                               (-CPB)  Size of the cloud-only prefetch buffer (in MB; 0 to disable).  Default value: 40.
clusterSize                          Optional<Integer>           --cluster-size                                        (-cluster)  The number of SNPs which make up a cluster. Must be at least 2  Default value: 3.
clusterWindowSize                    Optional<Integer>           --cluster-window-size                                 (-window)  The window size (in bases) in which to evaluate clustered SNPs  Default value: 0.
createOutputBamIndex                 Optional<Boolean>           --create-output-bam-index                             (-OBI)  If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.  Default value: true. Possible values: {true, false}
createOutputBamMd5                   Optional<Boolean>           --create-output-bam-md5                               (-OBM)  If true, create a MD5 digest for any BAM/SAM/CRAM file created  Default value: false. Possible values: {true, false}
createOutputVariantIndex             Optional<Boolean>           --create-output-variant-index                         (-OVI)  If true, create a VCF index when writing a coordinate-sorted VCF file.  Default value: true. Possible values: {true, false}
createOutputVariantMd5               Optional<Boolean>           --create-output-variant-md5                           (-OVM)  If true, create a a MD5 digest any VCF file created.  Default value: false. Possible values: {true, false}
disableBamIndexCaching               Optional<Boolean>           --disable-bam-index-caching                           (-DBIC)  If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.  Default value: false. Possible values: {true, false}
disableReadFilter                    Optional<String>            --disable-read-filter                                 (-DF)  Read filters to be disabled before analysis  This argument may be specified 0 or more times. Default value: null. Possible Values: {WellformedReadFilter}
disableSequenceDictionaryValidation  Optional<Boolean>           --disable-sequence-dictionary-validation              (-disable-sequence-dictionary-validation)  If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!  Default value: false. Possible values: {true, false}
excludeIntervals                     Optional<Boolean>           --exclude-intervals                                   (-XL) This argument may be specified 0 or more times. Default value: null.
filterExpression                     Optional<Array<String>>     --filter-expression                                   (-filter)  One or more expressions used with INFO fields to filter  This argument may be specified 0 or more times. Default value: null.
filterNotInMask                      Optional<Boolean>           --filter-not-in-mask                                  Filter records NOT in given input mask. Default value: false. Possible values: {true, false}
gatkConfigFile                       Optional<String>            --gatk-config-file                                    A configuration file to use with the GATK. Default value: null.
gcsMaxRetries                        Optional<Integer>           --gcs-max-retries                                     (-gcs-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20.
gcsProjectForRequesterPays           Optional<String>            --gcs-project-for-requester-pays                      Project to bill when accessing 'requester pays' buckets. If unset, these buckets cannot be accessed.  Default value: .
genotypeFilterExpression             Optional<String>            --genotype-filter-expression                          (-G-filter)  One or more expressions used with FORMAT (sample/genotype-level) fields to filter (see documentation guide for more info)  This argument may be specified 0 or more times. Default value: null.
genotypeFilterName                   Optional<String>            --genotype-filter-name                                (-G-filter-name)  Names to use for the list of sample/genotype filters (must be a 1-to-1 mapping); this name is put in the FILTER field for variants that get filtered  This argument may be specified 0 or more times. Default value: null.
help                                 Optional<Boolean>           --help                                                (-h) display the help message Default value: false. Possible values: {true, false}
inp                                  Optional<String>            --input                                               (-I) BAM/SAM/CRAM file containing reads This argument may be specified 0 or more times. Default value: null.
intervalExclusionPadding             Optional<Integer>           --interval-exclusion-padding                          (-ixp)  Amount of padding (in bp) to add to each interval you are excluding.  Default value: 0.
intervalMergingRule                  Optional<Boolean>           --interval-merging-rule                               (-imr)  Interval merging rule for abutting intervals  Default value: ALL. Possible values: {ALL, OVERLAPPING_ONLY}
intervalPadding                      Optional<Boolean>           --interval-padding                                    (-ip) Default value: 0.
intervalSetRule                      Optional<Boolean>           --interval-set-rule                                   (-isr)  Set merging approach to use for combining interval inputs  Default value: UNION. Possible values: {UNION, INTERSECTION}
intervals                            Optional<String>            --intervals                                           (-L) One or more genomic intervals over which to operate This argument may be specified 0 or more times. Default value: null.
invalidatePreviousFilters            Optional<Boolean>           --invalidate-previous-filters                         Remove previous filters applied to the VCF  Default value: false. Possible values: {true, false}
invertFilterExpression               Optional<Boolean>           --invert-filter-expression                            (-invfilter)  Invert the selection criteria for --filter-expression  Default value: false. Possible values: {true, false}
invertGenotypeFilterExpression       Optional<Boolean>           --invert-genotype-filter-expression                   (-invG-filter)  Invert the selection criteria for --genotype-filter-expression  Default value: false. Possible values: {true, false}
lenient                              Optional<Boolean>           --lenient                                             (-LE) Lenient processing of VCF files Default value: false. Possible values: {true, false}
mask                                 Optional<Boolean>           --mask                                                (-mask) Input mask Default value: null.
maskExtension                        Optional<Integer>           --mask-extension                                      How many bases beyond records from a provided 'mask' should variants be filtered Default value: 0.
maskName                             Optional<String>            --mask-name                                           The text to put in the FILTER field if a 'mask' is provided and overlaps with a variant call  Default value: Mask.
missingValuesEvaluateAsFailing       Optional<Boolean>           --missing-values-evaluate-as-failing                  When evaluating the JEXL expressions, missing values should be considered failing the expression  Default value: false. Possible values: {true, false}
quiet                                Optional<Boolean>           --QUIET                                               Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false}
readFilter                           Optional<String>            --read-filter                                         (-RF) Read filters to be applied before analysis This argument may be specified 0 or more times. Default value: null. Possible Values: {AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, SoftClippedReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter}
readIndex                            Optional<String>            --read-index                                          (-read-index)  Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically.  This argument may be specified 0 or more times. Default value: null.
readValidationStringency             Optional<Boolean>           --read-validation-stringency                          (-VS)  Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: SILENT. Possible values: {STRICT, LENIENT, SILENT}
reference                            Optional<FastaWithIndexes>  --reference                                           (-R) Reference sequence Default value: null.
secondsBetweenProgressUpdates        Optional<Double>            --seconds-between-progress-updates                    (-seconds-between-progress-updates)  Output traversal statistics every time this many seconds elapse  Default value: 10.0.
sequenceDictionary                   Optional<String>            --sequence-dictionary                                 (-sequence-dictionary)  Use the given sequence dictionary as the master/canonical sequence dictionary.  Must be a .dict file.  Default value: null.
setFilteredGenotypeToNoCall          Optional<Boolean>           --set-filtered-genotype-to-no-call                    Set filtered genotypes to no-call  Default value: false. Possible values: {true, false}
sitesOnlyVcfOutput                   Optional<Boolean>           --sites-only-vcf-output                               If true, don't emit genotype fields when writing vcf file output.  Default value: false. Possible values: {true, false}
tmpDir                               Optional<Boolean>           --tmp-dir                                             Temp directory to use. Default value: null.
useJdkDeflater                       Optional<Boolean>           --use-jdk-deflater                                    (-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false}
useJdkInflater                       Optional<Boolean>           --use-jdk-inflater                                    (-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false}
verbosity                            Optional<Boolean>           --verbosity                                           (-verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG}
version                              Optional<Boolean>           --version                                             display the version number for this tool Default value: false. Possible values: {true, false}
disableToolDefaultReadFilters        Optional<Boolean>           --disable-tool-default-read-filters                   (-disable-tool-default-read-filters)  Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)  Default value: false. Possible values: {true, false}
showhidden                           Optional<Boolean>           --showHidden                                          (-showHidden)  display hidden arguments  Default value: false. Possible values: {true, false}
ambigFilterBases                     Optional<Integer>           --ambig-filter-bases                                  Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction.  Default value: null.  Cannot be used in conjuction with argument(s) maxAmbiguousBaseFraction
ambigFilterFrac                      Optional<Double>            --ambig-filter-frac                                   Threshold fraction of ambiguous bases Default value: 0.05. Cannot be used in conjuction with argument(s) maxAmbiguousBases
maxFragmentLength                    Optional<Boolean>           --max-fragment-length                                 Default value: 1000000.
minFragmentLength                    Optional<Boolean>           --min-fragment-length                                 Default value: 0.
keepIntervals                        Optional<String>            --keep-intervals                                      One or more genomic intervals to keep This argument must be specified at least once. Required.
library                              Optional<String>            --library                                             (-library) Name of the library to keep This argument must be specified at least once. Required.
maximumMappingQuality                Optional<Integer>           --maximum-mapping-quality                             Maximum mapping quality to keep (inclusive)  Default value: null.
minimumMappingQuality                Optional<Integer>           --minimum-mapping-quality                             Minimum mapping quality to keep (inclusive)  Default value: 10.
dontRequireSoftClipsBothEnds         Optional<Boolean>           --dont-require-soft-clips-both-ends                   Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block  Default value: false. Possible values: {true, false}
filterTooShort                       Optional<Integer>           --filter-too-short                                    Minimum number of aligned bases Default value: 30.
platformFilterName                   Optional<Boolean>           --platform-filter-name                                This argument must be specified at least once. Required.
blackListedLanes                     Optional<String>            --black-listed-lanes                                  Platform unit (PU) to filter out This argument must be specified at least once. Required.
readGroupBlackList                   Optional<Boolean>           --read-group-black-list                               This argument must be specified at least once. Required.
keepReadGroup                        Optional<String>            --keep-read-group                                     The name of the read group to keep Required.
maxReadLength                        Optional<Integer>           --max-read-length                                     Keep only reads with length at most equal to the specified value Required.
minReadLength                        Optional<Integer>           --min-read-length                                     Keep only reads with length at least equal to the specified value Default value: 1.
readName                             Optional<String>            --read-name                                           Keep only reads with this read name Required.
keepReverseStrandOnly                Optional<Boolean>           --keep-reverse-strand-only                            Keep only reads on the reverse strand  Required. Possible values: {true, false}
sample                               Optional<String>            --sample                                              (-sample) The name of the sample(s) to keep, filtering out all others This argument must be specified at least once. Required.
invertSoftClipRatioFilter            Optional<Boolean>           --invert-soft-clip-ratio-filter                       Inverts the results from this filter, causing all variants that would pass to fail and visa-versa.  Default value: false. Possible values: {true, false}
softClippedLeadingTrailingRatio      Optional<Double>            --soft-clipped-leading-trailing-ratio                 Threshold ratio of soft clipped bases (leading / trailing the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumSoftClippedRatio
softClippedRatioThreshold            Optional<Double>            --soft-clipped-ratio-threshold                        Threshold ratio of soft clipped bases (anywhere in the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumLeadingTrailingSoftClippedRatio
===================================  ==========================  ========================================  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4Variantfiltration {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       String? outputFilename
       File variant
       File variant_tbi
       Boolean? addOutputSamProgramRecord
       Boolean? addOutputVcfCommandLine
       File? arguments_file
       Int? cloudIndexPrefetchBuffer
       Int? cloudPrefetchBuffer
       Int? clusterSize
       Int? clusterWindowSize
       Boolean? createOutputBamIndex
       Boolean? createOutputBamMd5
       Boolean? createOutputVariantIndex
       Boolean? createOutputVariantMd5
       Boolean? disableBamIndexCaching
       String? disableReadFilter
       Boolean? disableSequenceDictionaryValidation
       Boolean? excludeIntervals
       Array[String]? filterExpression
       Array[String?] filterName
       Boolean? filterNotInMask
       String? gatkConfigFile
       Int? gcsMaxRetries
       String? gcsProjectForRequesterPays
       String? genotypeFilterExpression
       String? genotypeFilterName
       Boolean? help
       String? inp
       Int? intervalExclusionPadding
       Boolean? intervalMergingRule
       Boolean? intervalPadding
       Boolean? intervalSetRule
       String? intervals
       Boolean? invalidatePreviousFilters
       Boolean? invertFilterExpression
       Boolean? invertGenotypeFilterExpression
       Boolean? lenient
       Boolean? mask
       Int? maskExtension
       String? maskName
       Boolean? missingValuesEvaluateAsFailing
       Boolean? quiet
       String? readFilter
       String? readIndex
       Boolean? readValidationStringency
       File? reference
       File? reference_fai
       File? reference_amb
       File? reference_ann
       File? reference_bwt
       File? reference_pac
       File? reference_sa
       File? reference_dict
       Float? secondsBetweenProgressUpdates
       String? sequenceDictionary
       Boolean? setFilteredGenotypeToNoCall
       Boolean? sitesOnlyVcfOutput
       Boolean? tmpDir
       Boolean? useJdkDeflater
       Boolean? useJdkInflater
       Boolean? verbosity
       Boolean? version
       Boolean? disableToolDefaultReadFilters
       Boolean? showhidden
       Int? ambigFilterBases
       Float? ambigFilterFrac
       Boolean? maxFragmentLength
       Boolean? minFragmentLength
       String? keepIntervals
       String? library
       Int? maximumMappingQuality
       Int? minimumMappingQuality
       Boolean? dontRequireSoftClipsBothEnds
       Int? filterTooShort
       Boolean? platformFilterName
       String? blackListedLanes
       Boolean? readGroupBlackList
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
       gatk VariantFiltration \
         --java-options '-Xmx~{((select_first([runtime_memory, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         --output '~{select_first([outputFilename, "~{basename(variant, ".vcf.gz")}.filtered.vcf"])}' \
         --variant '~{variant}' \
         ~{if (defined(addOutputSamProgramRecord) && select_first([addOutputSamProgramRecord])) then "--add-output-sam-program-record" else ""} \
         ~{if (defined(addOutputVcfCommandLine) && select_first([addOutputVcfCommandLine])) then "--add-output-vcf-command-line" else ""} \
         ~{if defined(arguments_file) then ("--arguments_file '" + arguments_file + "'") else ""} \
         ~{if defined(cloudIndexPrefetchBuffer) then ("--cloud-index-prefetch-buffer " + cloudIndexPrefetchBuffer) else ''} \
         ~{if defined(cloudPrefetchBuffer) then ("--cloud-prefetch-buffer " + cloudPrefetchBuffer) else ''} \
         ~{if defined(clusterSize) then ("--cluster-size " + clusterSize) else ''} \
         ~{if defined(clusterWindowSize) then ("--cluster-window-size " + clusterWindowSize) else ''} \
         ~{if (defined(createOutputBamIndex) && select_first([createOutputBamIndex])) then "--create-output-bam-index" else ""} \
         ~{if (defined(createOutputBamMd5) && select_first([createOutputBamMd5])) then "--create-output-bam-md5" else ""} \
         ~{if (defined(createOutputVariantIndex) && select_first([createOutputVariantIndex])) then "--create-output-variant-index" else ""} \
         ~{if (defined(createOutputVariantMd5) && select_first([createOutputVariantMd5])) then "--create-output-variant-md5" else ""} \
         ~{if (defined(disableBamIndexCaching) && select_first([disableBamIndexCaching])) then "--disable-bam-index-caching" else ""} \
         ~{if defined(disableReadFilter) then ("--disable-read-filter '" + disableReadFilter + "'") else ""} \
         ~{if (defined(disableSequenceDictionaryValidation) && select_first([disableSequenceDictionaryValidation])) then "--disable-sequence-dictionary-validation" else ""} \
         ~{if (defined(excludeIntervals) && select_first([excludeIntervals])) then "--exclude-intervals" else ""} \
         ~{if (defined(filterExpression) && length(select_first([filterExpression])) > 0) then "--filter-expression '" + sep("' --filter-expression '", select_first([filterExpression])) + "'" else ""} \
         ~{if length(filterName) > 0 then "--filter-name '" + sep("' --filter-name '", select_all(filterName)) + "'" else ""} \
         ~{if (defined(filterNotInMask) && select_first([filterNotInMask])) then "--filter-not-in-mask" else ""} \
         ~{if defined(gatkConfigFile) then ("--gatk-config-file '" + gatkConfigFile + "'") else ""} \
         ~{if defined(gcsMaxRetries) then ("--gcs-max-retries " + gcsMaxRetries) else ''} \
         ~{if defined(gcsProjectForRequesterPays) then ("--gcs-project-for-requester-pays '" + gcsProjectForRequesterPays + "'") else ""} \
         ~{if defined(genotypeFilterExpression) then ("--genotype-filter-expression '" + genotypeFilterExpression + "'") else ""} \
         ~{if defined(genotypeFilterName) then ("--genotype-filter-name '" + genotypeFilterName + "'") else ""} \
         ~{if (defined(help) && select_first([help])) then "--help" else ""} \
         ~{if defined(inp) then ("--input '" + inp + "'") else ""} \
         ~{if defined(intervalExclusionPadding) then ("--interval-exclusion-padding " + intervalExclusionPadding) else ''} \
         ~{if (defined(intervalMergingRule) && select_first([intervalMergingRule])) then "--interval-merging-rule" else ""} \
         ~{if (defined(intervalPadding) && select_first([intervalPadding])) then "--interval-padding" else ""} \
         ~{if (defined(intervalSetRule) && select_first([intervalSetRule])) then "--interval-set-rule" else ""} \
         ~{if defined(intervals) then ("--intervals '" + intervals + "'") else ""} \
         ~{if (defined(invalidatePreviousFilters) && select_first([invalidatePreviousFilters])) then "--invalidate-previous-filters" else ""} \
         ~{if (defined(invertFilterExpression) && select_first([invertFilterExpression])) then "--invert-filter-expression" else ""} \
         ~{if (defined(invertGenotypeFilterExpression) && select_first([invertGenotypeFilterExpression])) then "--invert-genotype-filter-expression" else ""} \
         ~{if (defined(lenient) && select_first([lenient])) then "--lenient" else ""} \
         ~{if (defined(mask) && select_first([mask])) then "--mask" else ""} \
         ~{if defined(maskExtension) then ("--mask-extension " + maskExtension) else ''} \
         ~{if defined(maskName) then ("--mask-name '" + maskName + "'") else ""} \
         ~{if (defined(missingValuesEvaluateAsFailing) && select_first([missingValuesEvaluateAsFailing])) then "--missing-values-evaluate-as-failing" else ""} \
         ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
         ~{if defined(readFilter) then ("--read-filter '" + readFilter + "'") else ""} \
         ~{if defined(readIndex) then ("--read-index '" + readIndex + "'") else ""} \
         ~{if (defined(readValidationStringency) && select_first([readValidationStringency])) then "--read-validation-stringency" else ""} \
         ~{if defined(reference) then ("--reference '" + reference + "'") else ""} \
         ~{if defined(secondsBetweenProgressUpdates) then ("--seconds-between-progress-updates " + secondsBetweenProgressUpdates) else ''} \
         ~{if defined(sequenceDictionary) then ("--sequence-dictionary '" + sequenceDictionary + "'") else ""} \
         ~{if (defined(setFilteredGenotypeToNoCall) && select_first([setFilteredGenotypeToNoCall])) then "--set-filtered-genotype-to-no-call" else ""} \
         ~{if (defined(sitesOnlyVcfOutput) && select_first([sitesOnlyVcfOutput])) then "--sites-only-vcf-output" else ""} \
         ~{if (defined(tmpDir) && select_first([tmpDir])) then "--tmp-dir" else ""} \
         ~{if (defined(useJdkDeflater) && select_first([useJdkDeflater])) then "--use-jdk-deflater" else ""} \
         ~{if (defined(useJdkInflater) && select_first([useJdkInflater])) then "--use-jdk-inflater" else ""} \
         ~{if (defined(verbosity) && select_first([verbosity])) then "--verbosity" else ""} \
         ~{if (defined(version) && select_first([version])) then "--version" else ""} \
         ~{if (defined(disableToolDefaultReadFilters) && select_first([disableToolDefaultReadFilters])) then "--disable-tool-default-read-filters" else ""} \
         ~{if (defined(showhidden) && select_first([showhidden])) then "--showHidden" else ""} \
         ~{if defined(ambigFilterBases) then ("--ambig-filter-bases " + ambigFilterBases) else ''} \
         ~{if defined(ambigFilterFrac) then ("--ambig-filter-frac " + ambigFilterFrac) else ''} \
         ~{if (defined(maxFragmentLength) && select_first([maxFragmentLength])) then "--max-fragment-length" else ""} \
         ~{if (defined(minFragmentLength) && select_first([minFragmentLength])) then "--min-fragment-length" else ""} \
         ~{if defined(keepIntervals) then ("--keep-intervals '" + keepIntervals + "'") else ""} \
         ~{if defined(library) then ("--library '" + library + "'") else ""} \
         ~{if defined(maximumMappingQuality) then ("--maximum-mapping-quality " + maximumMappingQuality) else ''} \
         ~{if defined(minimumMappingQuality) then ("--minimum-mapping-quality " + minimumMappingQuality) else ''} \
         ~{if (defined(dontRequireSoftClipsBothEnds) && select_first([dontRequireSoftClipsBothEnds])) then "--dont-require-soft-clips-both-ends" else ""} \
         ~{if defined(filterTooShort) then ("--filter-too-short " + filterTooShort) else ''} \
         ~{if (defined(platformFilterName) && select_first([platformFilterName])) then "--platform-filter-name" else ""} \
         ~{if defined(blackListedLanes) then ("--black-listed-lanes '" + blackListedLanes + "'") else ""} \
         ~{if (defined(readGroupBlackList) && select_first([readGroupBlackList])) then "--read-group-black-list" else ""} \
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
       docker: "broadinstitute/gatk:4.1.4.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "~{basename(variant, ".vcf.gz")}.filtered.vcf"])
       File out_tbi = select_first([outputFilename, "~{basename(variant, ".vcf.gz")}.filtered.vcf"]) + ".tbi"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: VariantFiltration'
   doc: |
     USAGE: VariantFiltration [arguments]
     Filter variant calls based on INFO and/or FORMAT annotations.
     Version:4.1.3.0

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
   - id: outputFilename
     label: outputFilename
     doc: (-O) File to which variants should be written Required.
     type:
     - string
     - 'null'
     default: generated.filtered.vcf
     inputBinding:
       prefix: --output
       valueFrom: $(inputs.variant.basename.replace(/.vcf.gz$/, "")).filtered.vcf
       separate: true
   - id: variant
     label: variant
     doc: (-V) A VCF file containing variants Required.
     type: File
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       prefix: --variant
       separate: true
   - id: addOutputSamProgramRecord
     label: addOutputSamProgramRecord
     doc: |-
       (-add-output-sam-program-record)  If true, adds a PG tag to created SAM/BAM/CRAM files.  Default value: true. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --add-output-sam-program-record
       separate: true
   - id: addOutputVcfCommandLine
     label: addOutputVcfCommandLine
     doc: |-
       (-add-output-vcf-command-line)  If true, adds a command line header line to created VCF files.  Default value: true. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --add-output-vcf-command-line
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
       (-CIPB)  Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.  Default value: -1. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --cloud-index-prefetch-buffer
       separate: true
   - id: cloudPrefetchBuffer
     label: cloudPrefetchBuffer
     doc: |-
       (-CPB)  Size of the cloud-only prefetch buffer (in MB; 0 to disable).  Default value: 40. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --cloud-prefetch-buffer
       separate: true
   - id: clusterSize
     label: clusterSize
     doc: |-
       (-cluster)  The number of SNPs which make up a cluster. Must be at least 2  Default value: 3. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --cluster-size
       separate: true
   - id: clusterWindowSize
     label: clusterWindowSize
     doc: |-
       (-window)  The window size (in bases) in which to evaluate clustered SNPs  Default value: 0. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --cluster-window-size
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
       separate: true
   - id: createOutputBamMd5
     label: createOutputBamMd5
     doc: |-
       (-OBM)  If true, create a MD5 digest for any BAM/SAM/CRAM file created  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --create-output-bam-md5
       separate: true
   - id: createOutputVariantIndex
     label: createOutputVariantIndex
     doc: |-
       (-OVI)  If true, create a VCF index when writing a coordinate-sorted VCF file.  Default value: true. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --create-output-variant-index
       separate: true
   - id: createOutputVariantMd5
     label: createOutputVariantMd5
     doc: |-
       (-OVM)  If true, create a a MD5 digest any VCF file created.  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --create-output-variant-md5
       separate: true
   - id: disableBamIndexCaching
     label: disableBamIndexCaching
     doc: |-
       (-DBIC)  If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --disable-bam-index-caching
       separate: true
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
       (-disable-sequence-dictionary-validation)  If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --disable-sequence-dictionary-validation
       separate: true
   - id: excludeIntervals
     label: excludeIntervals
     doc: '(-XL) This argument may be specified 0 or more times. Default value: null. '
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --exclude-intervals
       separate: true
   - id: filterExpression
     label: filterExpression
     doc: |-
       (-filter)  One or more expressions used with INFO fields to filter  This argument may be specified 0 or more times. Default value: null. 
     type:
     - type: array
       inputBinding:
         prefix: --filter-expression
         separate: true
       items: string
     - 'null'
     inputBinding: {}
   - id: filterName
     label: filterName
     doc: |-
       Names to use for the list of filters This argument may be specified 0 or more times. Default value: null. 
     type:
       type: array
       inputBinding:
         prefix: --filter-name
         separate: true
       items:
       - string
       - 'null'
     inputBinding: {}
   - id: filterNotInMask
     label: filterNotInMask
     doc: |-
       Filter records NOT in given input mask. Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --filter-not-in-mask
       separate: true
   - id: gatkConfigFile
     label: gatkConfigFile
     doc: 'A configuration file to use with the GATK. Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --gatk-config-file
       separate: true
   - id: gcsMaxRetries
     label: gcsMaxRetries
     doc: |-
       (-gcs-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --gcs-max-retries
       separate: true
   - id: gcsProjectForRequesterPays
     label: gcsProjectForRequesterPays
     doc: |2-
        Project to bill when accessing 'requester pays' buckets. If unset, these buckets cannot be accessed.  Default value: . 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --gcs-project-for-requester-pays
       separate: true
   - id: genotypeFilterExpression
     label: genotypeFilterExpression
     doc: |-
       (-G-filter)  One or more expressions used with FORMAT (sample/genotype-level) fields to filter (see documentation guide for more info)  This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --genotype-filter-expression
       separate: true
   - id: genotypeFilterName
     label: genotypeFilterName
     doc: |-
       (-G-filter-name)  Names to use for the list of sample/genotype filters (must be a 1-to-1 mapping); this name is put in the FILTER field for variants that get filtered  This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --genotype-filter-name
       separate: true
   - id: help
     label: help
     doc: |-
       (-h) display the help message Default value: false. Possible values: {true, false}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --help
       separate: true
   - id: inp
     label: inp
     doc: |-
       (-I) BAM/SAM/CRAM file containing reads This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --input
       separate: true
   - id: intervalExclusionPadding
     label: intervalExclusionPadding
     doc: |-
       (-ixp)  Amount of padding (in bp) to add to each interval you are excluding.  Default value: 0. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --interval-exclusion-padding
       separate: true
   - id: intervalMergingRule
     label: intervalMergingRule
     doc: |-
       (-imr)  Interval merging rule for abutting intervals  Default value: ALL. Possible values: {ALL, OVERLAPPING_ONLY} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --interval-merging-rule
       separate: true
   - id: intervalPadding
     label: intervalPadding
     doc: '(-ip) Default value: 0.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --interval-padding
       separate: true
   - id: intervalSetRule
     label: intervalSetRule
     doc: |-
       (-isr)  Set merging approach to use for combining interval inputs  Default value: UNION. Possible values: {UNION, INTERSECTION} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --interval-set-rule
       separate: true
   - id: intervals
     label: intervals
     doc: |-
       (-L) One or more genomic intervals over which to operate This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --intervals
       separate: true
   - id: invalidatePreviousFilters
     label: invalidatePreviousFilters
     doc: |2-
        Remove previous filters applied to the VCF  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --invalidate-previous-filters
       separate: true
   - id: invertFilterExpression
     label: invertFilterExpression
     doc: |-
       (-invfilter)  Invert the selection criteria for --filter-expression  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --invert-filter-expression
       separate: true
   - id: invertGenotypeFilterExpression
     label: invertGenotypeFilterExpression
     doc: |-
       (-invG-filter)  Invert the selection criteria for --genotype-filter-expression  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --invert-genotype-filter-expression
       separate: true
   - id: lenient
     label: lenient
     doc: |-
       (-LE) Lenient processing of VCF files Default value: false. Possible values: {true, false}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --lenient
       separate: true
   - id: mask
     label: mask
     doc: '(-mask) Input mask Default value: null.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --mask
       separate: true
   - id: maskExtension
     label: maskExtension
     doc: |-
       How many bases beyond records from a provided 'mask' should variants be filtered Default value: 0. 
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --mask-extension
       separate: true
   - id: maskName
     label: maskName
     doc: |-
       The text to put in the FILTER field if a 'mask' is provided and overlaps with a variant call  Default value: Mask. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --mask-name
       separate: true
   - id: missingValuesEvaluateAsFailing
     label: missingValuesEvaluateAsFailing
     doc: |2-
        When evaluating the JEXL expressions, missing values should be considered failing the expression  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --missing-values-evaluate-as-failing
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
       (-RF) Read filters to be applied before analysis This argument may be specified 0 or more times. Default value: null. Possible Values: {AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, SoftClippedReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter}
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read-filter
       separate: true
   - id: readIndex
     label: readIndex
     doc: |-
       (-read-index)  Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically.  This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read-index
       separate: true
   - id: readValidationStringency
     label: readValidationStringency
     doc: |-
       (-VS)  Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: SILENT. Possible values: {STRICT, LENIENT, SILENT} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --read-validation-stringency
       separate: true
   - id: reference
     label: reference
     doc: '(-R) Reference sequence Default value: null.'
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
   - id: secondsBetweenProgressUpdates
     label: secondsBetweenProgressUpdates
     doc: |-
       (-seconds-between-progress-updates)  Output traversal statistics every time this many seconds elapse  Default value: 10.0. 
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --seconds-between-progress-updates
       separate: true
   - id: sequenceDictionary
     label: sequenceDictionary
     doc: |-
       (-sequence-dictionary)  Use the given sequence dictionary as the master/canonical sequence dictionary.  Must be a .dict file.  Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sequence-dictionary
       separate: true
   - id: setFilteredGenotypeToNoCall
     label: setFilteredGenotypeToNoCall
     doc: |2-
        Set filtered genotypes to no-call  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --set-filtered-genotype-to-no-call
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
     - boolean
     - 'null'
     inputBinding:
       prefix: --tmp-dir
       separate: true
   - id: useJdkDeflater
     label: useJdkDeflater
     doc: |-
       (-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --use-jdk-deflater
       separate: true
   - id: useJdkInflater
     label: useJdkInflater
     doc: |-
       (-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --use-jdk-inflater
       separate: true
   - id: verbosity
     label: verbosity
     doc: |-
       (-verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --verbosity
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
       (-disable-tool-default-read-filters)  Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --disable-tool-default-read-filters
       separate: true
   - id: showhidden
     label: showhidden
     doc: |-
       (-showHidden)  display hidden arguments  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --showHidden
       separate: true
   - id: ambigFilterBases
     label: ambigFilterBases
     doc: |-
       Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction.  Default value: null.  Cannot be used in conjuction with argument(s) maxAmbiguousBaseFraction
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --ambig-filter-bases
       separate: true
   - id: ambigFilterFrac
     label: ambigFilterFrac
     doc: |-
       Threshold fraction of ambiguous bases Default value: 0.05. Cannot be used in conjuction with argument(s) maxAmbiguousBases
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --ambig-filter-frac
       separate: true
   - id: maxFragmentLength
     label: maxFragmentLength
     doc: 'Default value: 1000000.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --max-fragment-length
       separate: true
   - id: minFragmentLength
     label: minFragmentLength
     doc: 'Default value: 0.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --min-fragment-length
       separate: true
   - id: keepIntervals
     label: keepIntervals
     doc: |-
       One or more genomic intervals to keep This argument must be specified at least once. Required. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --keep-intervals
       separate: true
   - id: library
     label: library
     doc: |-
       (-library) Name of the library to keep This argument must be specified at least once. Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --library
       separate: true
   - id: maximumMappingQuality
     label: maximumMappingQuality
     doc: ' Maximum mapping quality to keep (inclusive)  Default value: null. '
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --maximum-mapping-quality
       separate: true
   - id: minimumMappingQuality
     label: minimumMappingQuality
     doc: ' Minimum mapping quality to keep (inclusive)  Default value: 10. '
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --minimum-mapping-quality
       separate: true
   - id: dontRequireSoftClipsBothEnds
     label: dontRequireSoftClipsBothEnds
     doc: |2-
        Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --dont-require-soft-clips-both-ends
       separate: true
   - id: filterTooShort
     label: filterTooShort
     doc: 'Minimum number of aligned bases Default value: 30.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --filter-too-short
       separate: true
   - id: platformFilterName
     label: platformFilterName
     doc: This argument must be specified at least once. Required.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --platform-filter-name
       separate: true
   - id: blackListedLanes
     label: blackListedLanes
     doc: |-
       Platform unit (PU) to filter out This argument must be specified at least once. Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --black-listed-lanes
       separate: true
   - id: readGroupBlackList
     label: readGroupBlackList
     doc: 'This argument must be specified at least once. Required. '
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --read-group-black-list
       separate: true
   - id: keepReadGroup
     label: keepReadGroup
     doc: The name of the read group to keep Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --keep-read-group
       separate: true
   - id: maxReadLength
     label: maxReadLength
     doc: Keep only reads with length at most equal to the specified value Required.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-read-length
       separate: true
   - id: minReadLength
     label: minReadLength
     doc: |-
       Keep only reads with length at least equal to the specified value Default value: 1.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-read-length
       separate: true
   - id: readName
     label: readName
     doc: Keep only reads with this read name Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read-name
       separate: true
   - id: keepReverseStrandOnly
     label: keepReverseStrandOnly
     doc: |2-
        Keep only reads on the reverse strand  Required. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --keep-reverse-strand-only
       separate: true
   - id: sample
     label: sample
     doc: |-
       (-sample) The name of the sample(s) to keep, filtering out all others This argument must be specified at least once. Required. 
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
       glob: $(inputs.variant.basename.replace(/.vcf.gz$/, "")).filtered.vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - VariantFiltration
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4Variantfiltration


