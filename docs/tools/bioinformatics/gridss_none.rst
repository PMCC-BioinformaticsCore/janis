:orphan:


Gridss
===============

Description
-------------

Tool identifier: ``gridss``

Tool path: ``janis_bioinformatics.tools.pappenfuss.gridss.versions import Gridss_2_2_3``

Version: None

Docker: ``gridss/gridss:v2.2.3``



Documentation
-------------

URL
******
*No URL to the documentation was provided*

Description
*********
*No documentation was provided: `contribute one <https://github.com/illusional>`_*

Outputs
-------
========  ======  ===============
name      type    documentation
========  ======  ===============
vcf       VCF
assembly  BAM
========  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

=========  =============  ===================  ==========  ===================================================================================================================================================================================================================================================================================================================================
name       type           prefix               position    documentation
=========  =============  ===================  ==========  ===================================================================================================================================================================================================================================================================================================================================
reference  FastaWithDict  REFERENCE_SEQUENCE=
bams       Array<BAM>     INPUT=                           (I=File Coordinate-sorted input BAM file. Default value: null. This option may be specified 0 or more times.
blacklist  bed            BLACKLIST=                       (BL=File) BED blacklist of regions to ignore. Assembly of regions such as high-coverage centromeric repeats is slow, and if such regions are to be filtered in downstream analysis anyway, blacklisting those region will improve runtime performance. For human WGS, the ENCODE DAC blacklist is recommended. Default value: null.
=========  =============  ===================  ==========  ===================================================================================================================================================================================================================================================================================================================================

Optional inputs
***************

=========================  ==================  =============================  ==========  =============================================================================================================================================================================================================================================================================================================
name                       type                prefix                         position    documentation
=========================  ==================  =============================  ==========  =============================================================================================================================================================================================================================================================================================================
outputFilename             Optional<Filename>  OUTPUT=                                    (O=) VCF structural variation calls. Required.
assemblyFilename           Optional<Filename>  ASSEMBLY=                                  Breakend assemblies which have undergone split read identification Required.
inputLabel                 Optional<String>    INPUT_LABEL=                               Input label. Variant calling evidence breakdowns are reported for each label. Default labels correspond to INPUT filenames. When specifying labels, labels must be provided for all input files. Default value: null. This option may be specified 0 or more times.
inputMaxFragmentSize       Optional<Integer>   INPUT_MAX_FRAGMENT_SIZE=                   Per input maximum concordant fragment size. Default value: null. This option may be specified 0 or more times.
inputMinFragmentSize       Optional<Integer>   INPUT_MIN_FRAGMENT_SIZE=                   Per input minimum concordant fragment size. Default value: null. This option may be specified 0 or more times.
readPairConcordantPercent  Optional<Float>     READ_PAIR_CONCORDANT_PERCENT=              Percent of read pairs considered concorant (0.0-1.0). If this is unset, the SAM proper pair flag is used to determine whether a read is discordantly aligned. Explicit fragment size specification overrides this setting. Default value: 0.995. This option can be set to 'null' to clear the default value.
configurationFile          Optional<File>      CONFIGURATION_FILE=                        (C=File) gridss configuration file containing overrides Default value: null.
workerThreads              Optional<Integer>   WORKER_THREADS=                            (THREADS=Integer  Number of worker threads to spawn. Defaults to number of cores available. Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread. Default value: 6. This option can be set to 'null' to clear the default value.
workingDir                 Optional<String>    WORKING_DIR=                               Directory to place intermediate results directories. Default location is the same directory as the associated input or output file. Default value: null.
ignoreDuplicates           Optional<Boolean>   IGNORE_DUPLICATES=                         Ignore reads marked as duplicates. Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false}
=========================  ==================  =============================  ==========  =============================================================================================================================================================================================================================================================================================================


Metadata
********

Author: **Unknown**


*Gridss was last updated on **Unknown***.
*This page was automatically generated on 2019-07-24*.
