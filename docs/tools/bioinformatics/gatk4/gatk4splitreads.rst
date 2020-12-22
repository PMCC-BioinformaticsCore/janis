:orphan:

GATK4: SplitReads
===================================

``Gatk4SplitReads`` · *1 contributor · 3 versions*

USAGE: SplitReads [arguments]
Outputs reads from a SAM/BAM/CRAM by read group, sample and library name
Version:4.1.3.0


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.splitreads.versions import Gatk4SplitReads_4_1_4

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4splitreads_step",
           Gatk4SplitReads_4_1_4(
               outputFilename=None,
               bam=None,
           )
       )
       wf.output("out", source=gatk4splitreads_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4SplitReads:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4SplitReads > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam




5. Run Gatk4SplitReads with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4SplitReads





Information
------------

:ID: ``Gatk4SplitReads``
:URL: *No URL to the documentation was provided*
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.4.0
:Authors: Michael Franklin
:Citations: None
:Created: 2019-09-16
:Updated: 2019-09-16


Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedBam  Bam
======  ==========  ===============


Additional configuration (inputs)
---------------------------------

===================================  ==========================  =======================================  ==========  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                                 type                        prefix                                     position  documentation
===================================  ==========================  =======================================  ==========  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
outputFilename                       String                      --output                                             The directory to output SAM/BAM/CRAM files. Default value: '.'
bam                                  IndexedBam                  --input                                           1  (-I:String) BAM/SAM/CRAM file containing reads  This argument must be specified at least once.
intervals                            Optional<bed>               --intervals                                          (-L:String) One or more genomic intervals over which to operate This argument may be specified 0 or more times. Default value: null.
javaOptions                          Optional<Array<String>>
compression_level                    Optional<Integer>                                                                Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
addOutputSamProgramRecord            Optional<Boolean>           -add-output-sam-program-record                       (--add-output-sam-program-record)  If true, adds a PG tag to created SAM/BAM/CRAM files.  Default value: true. Possible values: {true, false}
addOutputVcfCommandLine              Optional<Boolean>           -add-output-vcf-command-line                         (--add-output-vcf-command-line)  If true, adds a command line header line to created VCF files.  Default value: true. Possible values: {true, false}
arguments_file                       Optional<File>              --arguments_file:File                                read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null.
cloudIndexPrefetchBuffer             Optional<String>            --cloud-index-prefetch-buffer                        (-CIPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.  Default value: -1.
cloudPrefetchBuffer                  Optional<String>            --cloud-prefetch-buffer                              (-CPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable).  Default value: 40.
createOutputBamIndex                 Optional<String>            --create-output-bam-index                            (-OBI:Boolean)  If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.  Default value: true. Possible values: {true, false}
createOutputBamMd5                   Optional<String>            --create-output-bam-md5                              (-OBM:Boolean)  If true, create a MD5 digest for any BAM/SAM/CRAM file created  Default value: false. Possible values: {true, false}
createOutputVariantIndex             Optional<String>            --create-output-variant-index                        (-OVI:Boolean)  If true, create a VCF index when writing a coordinate-sorted VCF file.  Default value: true. Possible values: {true, false}
createOutputVariantMd5               Optional<String>            --create-output-variant-md5                          (-OVM:Boolean)  If true, create a a MD5 digest any VCF file created.  Default value: false. Possible values: {true, false}
disableBamIndexCaching               Optional<String>            --disable-bam-index-caching                          (-DBIC:Boolean)  If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.  Default value: false. Possible values: {true, false}
disableReadFilter                    Optional<String>            --disable-read-filter                                (-DF:String)  Read filters to be disabled before analysis  This argument may be specified 0 or more times. Default value: null. Possible Values: {WellformedReadFilter}
disableSequenceDictionaryValidation  Optional<Boolean>           -disable-sequence-dictionary-validation              (--disable-sequence-dictionary-validation)  If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!  Default value: false. Possible values: {true, false}
excludeIntervals                     Optional<String>            --exclude-intervals                                  (-XL:StringOne) This argument may be specified 0 or more times. Default value: null.
gatkConfigFile                       Optional<File>              --gatk-config-file                                   A configuration file to use with the GATK. Default value: null.
gcsRetries                           Optional<Integer>           -gcs-retries                                         (--gcs-max-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20.
gcsProjectForRequesterPays           Optional<String>            --gcs-project-for-requester-pays                     Project to bill when accessing requester pays  buckets. If unset, these buckets cannot be accessed.  Default value: .
intervalExclusionPadding             Optional<Integer>           --interval-exclusion-padding                         (-ixp:Integer)  Amount of padding (in bp) to add to each interval you are excluding.  Default value: 0.
imr                                  Optional<String>            -imr:IntervalMergingRule                             (--interval-merging-rule)  Interval merging rule for abutting intervals  Default value: ALL. Possible values: {ALL, OVERLAPPING_ONLY}
ip                                   Optional<Integer>           -ip                                                  (--interval-padding) Default value: 0.
isr                                  Optional<String>            -isr:IntervalSetRule                                 (--interval-set-rule)  Set merging approach to use for combining interval inputs  Default value: UNION. Possible values: {UNION, INTERSECTION}
le                                   Optional<Boolean>           --lenient                                            (-LE) Lenient processing of VCF files Default value: false. Possible values: {true, false}
quiet                                Optional<Boolean>           --QUIET                                              Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false}
readFilter                           Optional<String>            --read-filter                                        (-RF:String) Read filters to be applied before analysis This argument may be specified 0 or more times. Default value: null. Possible Values: {AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, SoftClippedReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter}
readIndex                            Optional<String>            -read-index                                          (--read-index)  Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically.  This argument may be specified 0 or more times. Default value: null.
readValidationStringency             Optional<String>            --read-validation-stringency                         (-VS:ValidationStringency)  Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: SITool returned: 0 LENT. Possible values: {STRICT, LENIENT, SILENT}
reference                            Optional<FastaWithIndexes>  --reference                                          (-R:String) Reference sequence Default value: null.
secondsBetweenProgressUpdates        Optional<Double>            -seconds-between-progress-updates                    (--seconds-between-progress-updates)  Output traversal statistics every time this many seconds elapse  Default value: 10.0.
sequenceDictionary                   Optional<String>            -sequence-dictionary                                 (--sequence-dictionary)  Use the given sequence dictionary as the master/canonical sequence dictionary.  Must be a .dict file.  Default value: null.
sitesOnlyVcfOutput                   Optional<Boolean>           --sites-only-vcf-output:Boolean                      If true, don't emit genotype fields when writing vcf file output.  Default value: false. Possible values: {true, false}
splitLibraryName                     Optional<String>            --split-library-name                                 (-LB)  Split file by library.  Default value: false. Possible values: {true, false}
rg                                   Optional<String>            --split-read-group                                   (-RG:BooleanSplit) Default value: false. Possible values: {true, false}
splitSample                          Optional<String>            --split-sample                                       (-SM:Boolean) Split file by sample. Default value: false. Possible values: {true, false}
tmpDir                               Optional<String>            --tmp-dir:GATKPathSpecifier                          Temp directory to use. Default value: null.
jdkDeflater                          Optional<Boolean>           -jdk-deflater                                        (--use-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false}
jdkInflater                          Optional<Boolean>           -jdk-inflater                                        (--use-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false}
verbosity                            Optional<String>            -verbosity:LogLevel                                  (--verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG}
disableToolDefaultReadFilters        Optional<Boolean>           -disable-tool-default-read-filters                   (--disable-tool-default-read-filters)  Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)  Default value: false. Possible values: {true, false}
ambigFilterBases                     Optional<Integer>           --ambig-filter-bases                                 Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction.  Default value: null.  Cannot be used in conjuction with argument(s) maxAmbiguousBaseFraction
ambigFilterFrac                      Optional<Double>            --ambig-filter-frac                                  Threshold fraction of ambiguous bases Default value: 0.05. Cannot be used in conjuction with argument(s) maxAmbiguousBases
maxFragmentLength                    Optional<Integer>           --max-fragment-length                                Default value: 1000000.
minFragmentLength                    Optional<Integer>           --min-fragment-length                                Default value: 0.
keepIntervals                        Optional<String>            --keep-intervals                                     Valid only if "IntervalOverlapReadFilter" is specified: One or more genomic intervals to keep This argument must be specified at least once. Required.
library                              Optional<String>            -library                                             (--library) Valid only if "LibraryReadFilter" is specified: Name of the library to keep This argument must be specified at least once. Required.
maximumMappingQuality                Optional<Integer>           --maximum-mapping-quality                            Maximum mapping quality to keep (inclusive)  Default value: null.
minimumMappingQuality                Optional<Integer>           --minimum-mapping-quality                            Minimum mapping quality to keep (inclusive)  Default value: 10.
dontRequireSoftClipsBothEnds         Optional<Boolean>           --dont-require-soft-clips-both-ends                  Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block  Default value: false. Possible values: {true, false}
filterTooShort                       Optional<Integer>           --filter-too-short                                   Minimum number of aligned bases Default value: 30.
platformFilterName                   Optional<String>            --platform-filter-name:String                        This argument must be specified at least once. Required.
blackListedLanes                     Optional<String>            --black-listed-lanes:String                          Platform unit (PU) to filter out This argument must be specified at least once. Required.
readGroupBlackList                   Optional<String>            --read-group-black-list:StringThe                    This argument must be specified at least once. Required.
keepReadGroup                        Optional<String>            --keep-read-group:String                             The name of the read group to keep Required.
maxReadLength                        Optional<Integer>           --max-read-length                                    Keep only reads with length at most equal to the specified value Required.
minReadLength                        Optional<Integer>           --min-read-length                                    Keep only reads with length at least equal to the specified value Default value: 1.
readName                             Optional<String>            --read-name:String                                   Keep only reads with this read name Required.
keepReverseStrandOnly                Optional<Boolean>           --keep-reverse-strand-only                           Keep only reads on the reverse strand  Required. Possible values: {true, false}
sample                               Optional<String>            -sample:String                                       (--sample) The name of the sample(s) to keep, filtering out all others This argument must be specified at least once. Required.
invertSoftClipRatioFilter            Optional<Boolean>           --invert-soft-clip-ratio-filter                      Inverts the results from this filter, causing all variants that would pass to fail and visa-versa.  Default value: false. Possible values: {true, false}
softClippedLeadingTrailingRatio      Optional<Double>            --soft-clipped-leading-trailing-ratio                Threshold ratio of soft clipped bases (leading / trailing the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumSoftClippedRatio
softClippedRatioThreshold            Optional<Double>            --soft-clipped-ratio-threshold                       Threshold ratio of soft clipped bases (anywhere in the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumLeadingTrailingSoftClippedRatio
===================================  ==========================  =======================================  ==========  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4SplitReads {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       String? outputFilename
       File bam
       File bam_bai
       File? intervals
       Array[String]? javaOptions
       Int? compression_level
       Boolean? addOutputSamProgramRecord
       Boolean? addOutputVcfCommandLine
       File? arguments_file
       String? cloudIndexPrefetchBuffer
       String? cloudPrefetchBuffer
       String? createOutputBamIndex
       String? createOutputBamMd5
       String? createOutputVariantIndex
       String? createOutputVariantMd5
       String? disableBamIndexCaching
       String? disableReadFilter
       Boolean? disableSequenceDictionaryValidation
       String? excludeIntervals
       File? gatkConfigFile
       Int? gcsRetries
       String? gcsProjectForRequesterPays
       Int? intervalExclusionPadding
       String? imr
       Int? ip
       String? isr
       Boolean? le
       Boolean? quiet
       String? readFilter
       String? readIndex
       String? readValidationStringency
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
       Boolean? sitesOnlyVcfOutput
       String? splitLibraryName
       String? rg
       String? splitSample
       String? tmpDir
       Boolean? jdkDeflater
       Boolean? jdkInflater
       String? verbosity
       Boolean? disableToolDefaultReadFilters
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
       cp -f '~{bam_bai}' $(echo '~{bam}' | sed 's/\.[^.]*$//').bai
       gatk SplitReads \
         --java-options '-Xmx~{((select_first([runtime_memory, 4, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         --output '~{select_first([outputFilename, "."])}' \
         ~{if defined(intervals) then ("--intervals '" + intervals + "'") else ""} \
         ~{if (defined(addOutputSamProgramRecord) && select_first([addOutputSamProgramRecord])) then "-add-output-sam-program-record" else ""} \
         ~{if (defined(addOutputVcfCommandLine) && select_first([addOutputVcfCommandLine])) then "-add-output-vcf-command-line" else ""} \
         ~{if defined(arguments_file) then ("--arguments_file:File '" + arguments_file + "'") else ""} \
         ~{if defined(cloudIndexPrefetchBuffer) then ("--cloud-index-prefetch-buffer '" + cloudIndexPrefetchBuffer + "'") else ""} \
         ~{if defined(cloudPrefetchBuffer) then ("--cloud-prefetch-buffer '" + cloudPrefetchBuffer + "'") else ""} \
         ~{if defined(createOutputBamIndex) then ("--create-output-bam-index '" + createOutputBamIndex + "'") else ""} \
         ~{if defined(createOutputBamMd5) then ("--create-output-bam-md5 '" + createOutputBamMd5 + "'") else ""} \
         ~{if defined(createOutputVariantIndex) then ("--create-output-variant-index '" + createOutputVariantIndex + "'") else ""} \
         ~{if defined(createOutputVariantMd5) then ("--create-output-variant-md5 '" + createOutputVariantMd5 + "'") else ""} \
         ~{if defined(disableBamIndexCaching) then ("--disable-bam-index-caching '" + disableBamIndexCaching + "'") else ""} \
         ~{if defined(disableReadFilter) then ("--disable-read-filter '" + disableReadFilter + "'") else ""} \
         ~{if (defined(disableSequenceDictionaryValidation) && select_first([disableSequenceDictionaryValidation])) then "-disable-sequence-dictionary-validation" else ""} \
         ~{if defined(excludeIntervals) then ("--exclude-intervals '" + excludeIntervals + "'") else ""} \
         ~{if defined(gatkConfigFile) then ("--gatk-config-file '" + gatkConfigFile + "'") else ""} \
         ~{if defined(gcsRetries) then ("-gcs-retries " + gcsRetries) else ''} \
         ~{if defined(gcsProjectForRequesterPays) then ("--gcs-project-for-requester-pays '" + gcsProjectForRequesterPays + "'") else ""} \
         ~{if defined(intervalExclusionPadding) then ("--interval-exclusion-padding " + intervalExclusionPadding) else ''} \
         ~{if defined(imr) then ("-imr:IntervalMergingRule '" + imr + "'") else ""} \
         ~{if defined(ip) then ("-ip " + ip) else ''} \
         ~{if defined(isr) then ("-isr:IntervalSetRule '" + isr + "'") else ""} \
         ~{if (defined(le) && select_first([le])) then "--lenient" else ""} \
         ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
         ~{if defined(readFilter) then ("--read-filter '" + readFilter + "'") else ""} \
         ~{if defined(readIndex) then ("-read-index '" + readIndex + "'") else ""} \
         ~{if defined(readValidationStringency) then ("--read-validation-stringency '" + readValidationStringency + "'") else ""} \
         ~{if defined(reference) then ("--reference '" + reference + "'") else ""} \
         ~{if defined(secondsBetweenProgressUpdates) then ("-seconds-between-progress-updates " + secondsBetweenProgressUpdates) else ''} \
         ~{if defined(sequenceDictionary) then ("-sequence-dictionary '" + sequenceDictionary + "'") else ""} \
         ~{if (defined(sitesOnlyVcfOutput) && select_first([sitesOnlyVcfOutput])) then "--sites-only-vcf-output:Boolean" else ""} \
         ~{if defined(splitLibraryName) then ("--split-library-name '" + splitLibraryName + "'") else ""} \
         ~{if defined(rg) then ("--split-read-group '" + rg + "'") else ""} \
         ~{if defined(splitSample) then ("--split-sample '" + splitSample + "'") else ""} \
         ~{if defined(tmpDir) then ("--tmp-dir:GATKPathSpecifier '" + tmpDir + "'") else ""} \
         ~{if (defined(jdkDeflater) && select_first([jdkDeflater])) then "-jdk-deflater" else ""} \
         ~{if (defined(jdkInflater) && select_first([jdkInflater])) then "-jdk-inflater" else ""} \
         ~{if defined(verbosity) then ("-verbosity:LogLevel '" + verbosity + "'") else ""} \
         ~{if (defined(disableToolDefaultReadFilters) && select_first([disableToolDefaultReadFilters])) then "-disable-tool-default-read-filters" else ""} \
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
         ~{if defined(platformFilterName) then ("--platform-filter-name:String '" + platformFilterName + "'") else ""} \
         ~{if defined(blackListedLanes) then ("--black-listed-lanes:String '" + blackListedLanes + "'") else ""} \
         ~{if defined(readGroupBlackList) then ("--read-group-black-list:StringThe '" + readGroupBlackList + "'") else ""} \
         ~{if defined(keepReadGroup) then ("--keep-read-group:String '" + keepReadGroup + "'") else ""} \
         ~{if defined(maxReadLength) then ("--max-read-length " + maxReadLength) else ''} \
         ~{if defined(minReadLength) then ("--min-read-length " + minReadLength) else ''} \
         ~{if defined(readName) then ("--read-name:String '" + readName + "'") else ""} \
         ~{if (defined(keepReverseStrandOnly) && select_first([keepReverseStrandOnly])) then "--keep-reverse-strand-only" else ""} \
         ~{if defined(sample) then ("-sample:String '" + sample + "'") else ""} \
         ~{if (defined(invertSoftClipRatioFilter) && select_first([invertSoftClipRatioFilter])) then "--invert-soft-clip-ratio-filter" else ""} \
         ~{if defined(softClippedLeadingTrailingRatio) then ("--soft-clipped-leading-trailing-ratio " + softClippedLeadingTrailingRatio) else ''} \
         ~{if defined(softClippedRatioThreshold) then ("--soft-clipped-ratio-threshold " + softClippedRatioThreshold) else ''} \
         --input '~{bam}'
       if [ -f $(echo '~{basename(bam)}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{basename(bam)}' | sed 's/\.[^.]*$//').bai $(echo '~{basename(bam)}' ).bai; fi
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.4.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4, 4])}G"
       preemptible: 2
     }
     output {
       File out = basename(bam)
       File out_bai = basename(bam) + ".bai"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: SplitReads'
   doc: |-
     USAGE: SplitReads [arguments]
     Outputs reads from a SAM/BAM/CRAM by read group, sample and library name
     Version:4.1.3.0

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.4.0

   inputs:
   - id: outputFilename
     label: outputFilename
     doc: "The directory to output SAM/BAM/CRAM files. Default value: '.' "
     type: string
     default: .
     inputBinding:
       prefix: --output
   - id: bam
     label: bam
     doc: |-
       (-I:String) BAM/SAM/CRAM file containing reads  This argument must be specified at least once.
     type: File
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
                           location: resolveSecondary(self.location, "^.bai"),
                           basename: resolveSecondary(self.basename, ".bai"),
                           class: "File",
                       }
               ];

       }
     inputBinding:
       prefix: --input
       position: 1
   - id: intervals
     label: intervals
     doc: |-
       (-L:String) One or more genomic intervals over which to operate This argument may be specified 0 or more times. Default value: null. 
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --intervals
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
   - id: arguments_file
     label: arguments_file
     doc: |-
       read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null. 
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --arguments_file:File
   - id: cloudIndexPrefetchBuffer
     label: cloudIndexPrefetchBuffer
     doc: |-
       (-CIPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.  Default value: -1. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --cloud-index-prefetch-buffer
   - id: cloudPrefetchBuffer
     label: cloudPrefetchBuffer
     doc: |-
       (-CPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable).  Default value: 40. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --cloud-prefetch-buffer
   - id: createOutputBamIndex
     label: createOutputBamIndex
     doc: |-
       (-OBI:Boolean)  If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.  Default value: true. Possible values: {true, false} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --create-output-bam-index
   - id: createOutputBamMd5
     label: createOutputBamMd5
     doc: |-
       (-OBM:Boolean)  If true, create a MD5 digest for any BAM/SAM/CRAM file created  Default value: false. Possible values: {true, false} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --create-output-bam-md5
   - id: createOutputVariantIndex
     label: createOutputVariantIndex
     doc: |-
       (-OVI:Boolean)  If true, create a VCF index when writing a coordinate-sorted VCF file.  Default value: true. Possible values: {true, false} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --create-output-variant-index
   - id: createOutputVariantMd5
     label: createOutputVariantMd5
     doc: |-
       (-OVM:Boolean)  If true, create a a MD5 digest any VCF file created.  Default value: false. Possible values: {true, false} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --create-output-variant-md5
   - id: disableBamIndexCaching
     label: disableBamIndexCaching
     doc: |-
       (-DBIC:Boolean)  If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.  Default value: false. Possible values: {true, false} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --disable-bam-index-caching
   - id: disableReadFilter
     label: disableReadFilter
     doc: |-
       (-DF:String)  Read filters to be disabled before analysis  This argument may be specified 0 or more times. Default value: null. Possible Values: {WellformedReadFilter}
     type:
     - string
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
   - id: excludeIntervals
     label: excludeIntervals
     doc: |-
       (-XL:StringOne) This argument may be specified 0 or more times. Default value: null. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --exclude-intervals
   - id: gatkConfigFile
     label: gatkConfigFile
     doc: 'A configuration file to use with the GATK. Default value: null.'
     type:
     - File
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
        Project to bill when accessing requester pays  buckets. If unset, these buckets cannot be accessed.  Default value: . 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --gcs-project-for-requester-pays
   - id: intervalExclusionPadding
     label: intervalExclusionPadding
     doc: |-
       (-ixp:Integer)  Amount of padding (in bp) to add to each interval you are excluding.  Default value: 0. 
     type:
     - int
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
       prefix: -imr:IntervalMergingRule
   - id: ip
     label: ip
     doc: '(--interval-padding) Default value: 0.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -ip
   - id: isr
     label: isr
     doc: |-
       (--interval-set-rule)  Set merging approach to use for combining interval inputs  Default value: UNION. Possible values: {UNION, INTERSECTION} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -isr:IntervalSetRule
   - id: le
     label: le
     doc: |-
       (-LE) Lenient processing of VCF files Default value: false. Possible values: {true, false}
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --lenient
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
       (-RF:String) Read filters to be applied before analysis This argument may be specified 0 or more times. Default value: null. Possible Values: {AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, SoftClippedReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter}
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
       (-VS:ValidationStringency)  Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: SITool returned: 0 LENT. Possible values: {STRICT, LENIENT, SILENT} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read-validation-stringency
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
       prefix: --sites-only-vcf-output:Boolean
   - id: splitLibraryName
     label: splitLibraryName
     doc: |-
       (-LB)  Split file by library.  Default value: false. Possible values: {true, false} 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --split-library-name
   - id: rg
     label: rg
     doc: '(-RG:BooleanSplit) Default value: false. Possible values: {true, false}'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --split-read-group
   - id: splitSample
     label: splitSample
     doc: |-
       (-SM:Boolean) Split file by sample. Default value: false. Possible values: {true, false}
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --split-sample
   - id: tmpDir
     label: tmpDir
     doc: 'Temp directory to use. Default value: null.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --tmp-dir:GATKPathSpecifier
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
       prefix: -verbosity:LogLevel
   - id: disableToolDefaultReadFilters
     label: disableToolDefaultReadFilters
     doc: |-
       (--disable-tool-default-read-filters)  Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -disable-tool-default-read-filters
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
       Valid only if "IntervalOverlapReadFilter" is specified: One or more genomic intervals to keep This argument must be specified at least once. Required. 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --keep-intervals
   - id: library
     label: library
     doc: |-
       (--library) Valid only if "LibraryReadFilter" is specified: Name of the library to keep This argument must be specified at least once. Required.
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
     doc: ' Minimum mapping quality to keep (inclusive)  Default value: 10. '
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
       prefix: --platform-filter-name:String
   - id: blackListedLanes
     label: blackListedLanes
     doc: |-
       Platform unit (PU) to filter out This argument must be specified at least once. Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --black-listed-lanes:String
   - id: readGroupBlackList
     label: readGroupBlackList
     doc: 'This argument must be specified at least once. Required. '
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --read-group-black-list:StringThe
   - id: keepReadGroup
     label: keepReadGroup
     doc: The name of the read group to keep Required.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --keep-read-group:String
   - id: maxReadLength
     label: maxReadLength
     doc: Keep only reads with length at most equal to the specified value Required.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-read-length
   - id: minReadLength
     label: minReadLength
     doc: |-
       Keep only reads with length at least equal to the specified value Default value: 1.
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
       prefix: --read-name:String
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
       prefix: -sample:String
   - id: invertSoftClipRatioFilter
     label: invertSoftClipRatioFilter
     doc: |2-
        Inverts the results from this filter, causing all variants that would pass to fail and visa-versa.  Default value: false. Possible values: {true, false} 
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --invert-soft-clip-ratio-filter
   - id: softClippedLeadingTrailingRatio
     label: softClippedLeadingTrailingRatio
     doc: |2-
        Threshold ratio of soft clipped bases (leading / trailing the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumSoftClippedRatio
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --soft-clipped-leading-trailing-ratio
   - id: softClippedRatioThreshold
     label: softClippedRatioThreshold
     doc: |2-
        Threshold ratio of soft clipped bases (anywhere in the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumLeadingTrailingSoftClippedRatio
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --soft-clipped-ratio-threshold

   outputs:
   - id: out
     label: out
     doc: Bam
     type: File
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
       glob: $(inputs.bam.basename)
       outputEval: $(inputs.bam.basename.basename)
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - SplitReads
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 4, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4SplitReads


