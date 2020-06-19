:orphan:

GATK4: SelectVariants
===========================================

*0 contributors · 4 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


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
:Authors: 
:Citations: None
:Created: None
:Updated: None



Outputs
-----------

======  ====================  ===============
name    type                  documentation
======  ====================  ===============
out     CompressedIndexedVCF
======  ====================  ===============



Additional configuration (inputs)
---------------------------------

===================================  ==============================  =======================================  ==========  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                                 type                            prefix                                   position    documentation
===================================  ==============================  =======================================  ==========  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
outputFilename                       Optional<Filename>              -O                                                   (--output) Required.
variants                             Optional<CompressedIndexedVCF>  -V                                                   (--variant) A VCF file containing variants Required.
addOutputSamProgramRecord            Optional<Boolean>               -add-output-sam-program-record:Boolean               (--add-output-sam-program-record)  If true, adds a PG tag to created SAM/BAM/CRAM files.  Default value: true. Possible values: {true, false}
addOutputVcfCommandLine              Optional<Boolean>               -add-output-vcf-command-line                         (--add-output-vcf-command-line)  If true, adds a command line header line to created VCF files.  Default value: true. Possible values: {true, false}
arguments_file                       Optional<File>                  --arguments_file                                     read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null.
cloudIndexPrefetchBuffer             Optional<Integer>               --cloud-index-prefetch-buffer                        (-CIPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.  Default value: -1.
cloudPrefetchBuffer                  Optional<Integer>               --cloud-prefetch-buffer                              (-CPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable).  Default value: 40.
conc                                 Optional<String>                -conc                                                (--concordance)  Output variants also called in this comparison track  Default value: null.
createOutputBamIndex                 Optional<Boolean>               --create-output-bam-index                            (-OBI)  If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.  Default value: true. Possible values: {true, false}
createOutputBamMd5                   Optional<Boolean>               --create-output-bam-md5                              (-OBM)  If true, create a MD5 digest for any BAM/SAM/CRAM file created  Default value: false. Possible values: {true, false}
createOutputvariantIndex             Optional<Boolean>               --create-output-variant-index                        (-OVI)  If true, create a VCF index when writing a coordinate-sorted VCF file.  Default value: true. Possible values: {true, false}
createOutputvariantMd5               Optional<Boolean>               --create-output-variant-md5                          (-OVM)  If true, create a a MD5 digest any VCF file created.  Default value: false. Possible values: {true, false}
disableBamIndexCaching               Optional<Boolean>               --disable-bam-index-caching                          (-DBIC:Boolean)  If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.  Default value: false. Possible values: {true, false}
disableReadFilter                    Optional<String>                --disable-read-filter                                (-DF)  Read filters to be disabled before analysis  This argument may be specified 0 or more times. Default value: null. Possible Values: {WellformedReadFilter}
disableSequenceDictionaryValidation  Optional<Boolean>               -disable-sequence-dictionary-validation              (--disable-sequence-dictionary-validation)  If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!  Default value: false. Possible values: {true, false}
disc                                 Optional<String>                -disc                                                (--discordance)  Output variants not called in this comparison track  Default value: null.
dropGenotypeAnnotation               Optional<String>                --drop-genotype-annotation                           (-DGA:String)  Genotype annotations to drop from output vcf.  Annotations to be dropped are specified by their key.  This argument may be specified 0 or more times. Default value: null.
dropInfoAnnotation                   Optional<String>                --drop-info-annotation                               (-DA:String)  Info annotations to drop from output vcf.  Annotations to be dropped are specified by their key.  This argument may be specified 0 or more times. Default value: null.
excludeFiltered                      Optional<Boolean>               --exclude-filtered                                   Don't include filtered sites Default value: false. Possible values: {true, false}
xlIds                                Optional<String>                -xl-ids                                              (--exclude-ids) List of variant rsIDs to exclude This argument may be specified 0 or more times. Default value: null.
excludeIntervals                     Optional<String>                --exclude-intervals                                  (-XL) This argument may be specified 0 or more times. Default value: null.
excludeNonvariants                   Optional<String>                --exclude-non-variants                               Default value: false. Possible values: {true, false}
excludeSampleExpressions             Optional<String>                --exclude-sample-expressions                         (-xl-se:String)  List of sample expressions to exclude  This argument may be specified 0 or more times. Default value: null.
excludeSampleName                    Optional<String>                --exclude-sample-name                                (-xl-sn:String)  Exclude genotypes from this sample  This argument may be specified 0 or more times. Default value: null.
gatkConfigFile                       Optional<File>                  --gatk-config-file                                   A configuration file to use with the GATK. Default value: null.
gcsRetries                           Optional<Integer>               -gcs-retries                                         (--gcs-max-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20.
gcsProjectForRequesterPays           Optional<String>                --gcs-project-for-requester-pays                     Project to bill when accessing requester pays buckets. If unset, these buckets cannot be accessed.  Default value: .
help                                 Optional<Boolean>               -h                                                   (--help) display the help message Default value: false. Possible values: {true, false}
bam                                  Optional<IndexedBam>            -I                                                   (--input) BAM/SAM/CRAM file containing reads This argument may be specified 0 or more times. Default value: null.
intervalExclusionPadding             Optional<Integer>               --interval-exclusion-padding                         (-ixp:Integer)  Amount of padding (in bp) to add to each interval you are excluding.  Default value: 0.
imr                                  Optional<String>                -imr                                                 (--interval-merging-rule)  Interval merging rule for abutting intervals  Default value: ALL. Possible values: {ALL, OVERLAPPING_ONLY}
ip                                   Optional<Integer>               -ip                                                  (--interval-padding) Default value: 0.
isr                                  Optional<String>                -isr                                                 (--interval-set-rule)  Set merging approach to use for combining interval inputs  Default value: UNION. Possible values: {UNION, INTERSECTION}
intervals                            Optional<String>                --intervals                                          (-L:String) One or more genomic intervals over which to operate This argument may be specified 0 or more times. Default value: null.
invertMendelianViolation             Optional<Boolean>               --invert-mendelian-violation                         Output non-mendelian violation sites only  Default value: false. Possible values: {true, false}
invertSelect                         Optional<Boolean>               -invert-select                                       (--invertSelect)  Invert the selection criteria for -select  Default value: false. Possible values: {true, false}
ids                                  Optional<String>                -ids                                                 (--keep-ids) List of variant rsIDs to select This argument may be specified 0 or more times. Default value: null.
keepOriginalAc                       Optional<Boolean>               --keep-original-ac                                   Store the original AC, AF, and AN values after subsetting Default value: false. Possible values: {true, false}
keepOriginalDp                       Optional<Boolean>               --keep-original-dp                                   Store the original DP value after subsetting Default value: false. Possible values: {true, false}
le                                   Optional<Boolean>               -LE                                                  (--lenient) Lenient processing of VCF files Default value: false. Possible values: {true, false}
maxFilteredGenotypes                 Optional<Integer>               --max-filtered-genotypes                             Maximum number of samples filtered at the genotype level  Default value: 2147483647.
maxFractionFilteredGenotypes         Optional<Double>                --max-fraction-filtered-genotypes                    Maximum fraction of samples filtered at the genotype level  Default value: 1.0.
maxIndelSize                         Optional<Integer>               --max-indel-size                                     Maximum size of indels to include Default value: 2147483647.
maxNocallFraction                    Optional<Double>                --max-nocall-fraction                                Maximum fraction of samples with no-call genotypes Default value: 1.0.
maxNocallNumber                      Optional<Integer>               --max-nocall-number                                  Maximum number of samples with no-call genotypes Default value: 2147483647.
mendelianViolation                   Optional<Boolean>               --mendelian-violation                                Default value: false. Possible values: {true, false}
mendelianViolationQualThreshold      Optional<Double>                --mendelian-violation-qual-threshold                 Minimum GQ score for each trio member to accept a site as a violation  Default value: 0.0.
minFilteredGenotypes                 Optional<Integer>               --min-filtered-genotypes                             Minimum number of samples filtered at the genotype level  Default value: 0.
minFractionFilteredGenotypes         Optional<Double>                --min-fraction-filtered-genotypes                    Maximum fraction of samples filtered at the genotype level  Default value: 0.0.
minIndelSize                         Optional<Integer>               --min-indel-size                                     Minimum size of indels to include Default value: 0.
pedigree                             Optional<File>                  --pedigree                                           (-ped:File) Pedigree file Default value: null.
preserveAlleles                      Optional<Boolean>               --preserve-alleles                                   Preserve original alleles, do not trim Default value: false. Possible values: {true, false}
quiet                                Optional<Boolean>               --QUIET                                              Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false}
readFilter                           Optional<String>                --read-filter                                        (-RF:String) Read filters to be applied before analysis This argument may be specified 0 or more times. Default value: null. Possible Values: {AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, SoftClippedReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter}
readIndex                            Optional<File>                  -read-index                                          (--read-index)  Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically.  This argument may be specified 0 or more times. Default value: null.
readValidationStringency             Optional<String>                --read-validation-stringency                         (-VS:ValidationStringency)  Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: SILENT. Possible values: {STRICT, LENIENT, SILENT}
reference                            Optional<FastaWithIndexes>      --reference                                          (-R:String) Reference sequence Default value: null.
removeFractionGenotypes              Optional<Double>                --remove-fraction-genotypes                          Select a fraction of genotypes at random from the input and sets them to no-call  Default value: 0.0.
removeUnusedAlternates               Optional<Boolean>               --remove-unused-alternates                           Remove alternate alleles not present in any genotypes  Default value: false. Possible values: {true, false}
restrictAllelesTo                    Optional<String>                --restrict-alleles-to                                Select only variants of a particular allelicity  Default value: ALL. Possible values: {ALL, BIALLELIC, MULTIALLELIC}
sampleExpressions                    Optional<String>                --sample-expressions                                 (-se:String)  Regular expression to select multiple samples  This argument may be specified 0 or more times. Default value: null.
sampleName                           Optional<String>                --sample-name                                        (-sn:String) Include genotypes from this sample This argument may be specified 0 or more times. Default value: null.
secondsBetweenProgressUpdates        Optional<Double>                -seconds-between-progress-updates                    (--seconds-between-progress-updates)  Output traversal statistics every time this many seconds elapse  Default value: 10.0.
selectRandomFraction                 Optional<String>                --select-random-fraction                             (-fraction:Double)  Select a fraction of variants at random from the input  Default value: 0.0.
selectTypeToExclude                  Optional<String>                --select-type-to-exclude                             (-xl-select-type:Type)  Do not select certain type of variants from the input file  This argument may be specified 0 or more times. Default value: null. Possible values: {NO_VARIATION, SNP, MNP, INDEL, SYMBOLIC, MIXED}
selectTypeToInclude                  Optional<String>                --select-type-to-include                             (-select-type:Type)  Select only a certain type of variants from the input file  This argument may be specified 0 or more times. Default value: null. Possible values: {NO_VARIATION, SNP, MNP, INDEL, SYMBOLIC, MIXED}
selectexpressions                    Optional<String>                --selectExpressions                                  (-select:String)  One or more criteria to use when selecting the data  This argument may be specified 0 or more times. Default value: null.
sequenceDictionary                   Optional<File>                  -sequence-dictionary                                 (--sequence-dictionary)  Use the given sequence dictionary as the master/canonical sequence dictionary.  Must be a .dict file.  Default value: null.
setFilteredGtToNocall                Optional<Boolean>               --set-filtered-gt-to-nocall                          Set filtered genotypes to no-call  Default value: false. Possible values: {true, false}
sitesOnlyVcfOutput                   Optional<Boolean>               --sites-only-vcf-output                              If true, don't emit genotype fields when writing vcf file output.  Default value: false. Possible values: {true, false}
tmpDir                               Optional<Filename>              --tmp-dir                                            Temp directory to use. Default value: null.
jdkDeflater                          Optional<Boolean>               -jdk-deflater                                        (--use-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false}
jdkInflater                          Optional<Boolean>               -jdk-inflater                                        (--use-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false}
verbosity                            Optional<String>                -verbosity                                           (--verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG}
version                              Optional<Boolean>               --version                                            display the version number for this tool Default value: false. Possible values: {true, false}
disableToolDefaultReadFilters        Optional<Boolean>               -disable-tool-default-read-filters                   (--disable-tool-default-read-filters)  Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)  Default value: false. Possible values: {true, false}
showhidden                           Optional<Boolean>               -showHidden                                          (--showHidden)  display hidden arguments  Default value: false. Possible values: {true, false}
ambigFilterBases                     Optional<Integer>               --ambig-filter-bases                                 Valid only if 'AmbiguousBaseReadFilter' is specified: Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction.  Default value: null.  Cannot be used in conjuction with argument(s) maxAmbiguousBaseFraction
ambigFilterFrac                      Optional<Double>                --ambig-filter-frac                                  Valid only if 'AmbiguousBaseReadFilter' is specified: Threshold fraction of ambiguous bases Default value: 0.05. Cannot be used in conjuction with argument(s) maxAmbiguousBases
maxFragmentLength                    Optional<Integer>               --max-fragment-length                                Valid only if 'FragmentLengthReadFilter' is specified: Maximum length of fragment (insert size) Default value: 1000000.
minFragmentLength                    Optional<Integer>               --min-fragment-length                                Valid only if 'FragmentLengthReadFilter' is specified: Minimum length of fragment (insert size) Default value: 0.
keepIntervals                        Optional<String>                --keep-intervals                                     Valid only if 'IntervalOverlapReadFilter' is specified: One or more genomic intervals to keep This argument must be specified at least once. Required.
library                              Optional<String>                -library                                             Valid only if 'LibraryReadFilter' is specified: (--library) Name of the library to keep This argument must be specified at least once. Required.
maximumMappingQuality                Optional<Integer>               --maximum-mapping-quality                            Valid only if 'MappingQualityReadFilter' is specified: Maximum mapping quality to keep (inclusive)  Default value: null.
minimumMappingQuality                Optional<Integer>               --minimum-mapping-quality                            Valid only if 'MappingQualityReadFilter' is specified: Minimum mapping quality to keep (inclusive)  Default value: 10.
dontRequireSoftClipsBothEnds         Optional<Boolean>               --dont-require-soft-clips-both-ends                  Valid only if 'OverclippedReadFilter' is specified: Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block  Default value: false. Possible values: {true, false}
filterTooShort                       Optional<Integer>               --filter-too-short                                   Valid only if 'OverclippedReadFilter' is specified: Minimum number of aligned bases Default value: 30.
platformFilterName                   Optional<String>                --platform-filter-name                               Valid only if 'PlatformReadFilter' is specified: This argument must be specified at least once. Required.
blackListedLanes                     Optional<String>                --black-listed-lanes                                 Valid only if 'PlatformUnitReadFilter' is specified: Platform unit (PU) to filter out This argument must be specified at least once. Required.
readGroupBlackList                   Optional<String>                --read-group-black-list                              Valid only if 'ReadGroupBlackListReadFilter' is specified: The name of the read group to filter out. This argument must be specified at least once. Required.
keepReadGroup                        Optional<String>                --keep-read-group                                    Valid only if 'ReadGroupReadFilter' is specified: The name of the read group to keep Required.
maxReadLength                        Optional<Integer>               --max-read-length                                    Valid only if 'ReadLengthReadFilter' is specified: Keep only reads with length at most equal to the specified value Required.
minReadLength                        Optional<Integer>               --min-read-length                                    Valid only if 'ReadLengthReadFilter' is specified: Keep only reads with length at least equal to the specified value Default value: 1.
readName                             Optional<String>                --read-name                                          Valid only if 'ReadNameReadFilter' is specified: Keep only reads with this read name Required.
keepReverseStrandOnly                Optional<Boolean>               --keep-reverse-strand-only                           Valid only if 'ReadStrandFilter' is specified: Keep only reads on the reverse strand  Required. Possible values: {true, false}
sample                               Optional<String>                --sample                                             Valid only if 'SampleReadFilter' is specified: The name of the sample(s) to keep, filtering out all others This argument must be specified at least once. Required.
invertSoftClipRatioFilter            Optional<Boolean>               --invert-soft-clip-ratio-filter                      Inverts the results from this filter, causing all variants that would pass to fail and visa-versa.  Default value: false. Possible values: {true, false}
softClippedLeadingTrailingRatio      Optional<Double>                --soft-clipped-leading-trailing-ratio                Threshold ratio of soft clipped bases (leading / trailing the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumSoftClippedRatio
softClippedRatioThreshold            Optional<Double>                --soft-clipped-ratio-threshold                       Threshold ratio of soft clipped bases (anywhere in the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumLeadingTrailingSoftClippedRatio
===================================  ==============================  =======================================  ==========  ======================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
