
.. include:: gatk4sortsam_4.0.12.0

GATK4: SortSAM
=============================

Description
-------------

Tool identifier: ``gatk4sortsam``

Tool path: ``janis_bioinformatics.tools.gatk4.sortsam.sortsam_4_0 import Gatk4SortSam_4_0``

Version: 4.0.12.0

Container: ``broadinstitute/gatk:4.0.12.0``



Documentation
-------------

URL
******
`https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.3/org_broadinstitute_hellbender_tools_picard_sam_SortSam.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.3/org_broadinstitute_hellbender_tools_picard_sam_SortSam.php>`_

Tool documentation
******************
Sorts a SAM/BAM/CRAM file.

Outputs
-------
======  =======  ===============
name    type     documentation
======  =======  ===============
out     BamPair
======  =======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

=========  ======  ========  ==========  ==============================================================================================================================================================
name       type    prefix      position  documentation
=========  ======  ========  ==========  ==============================================================================================================================================================
bam        BAM     -I                10  The SAM/BAM/CRAM file to sort.
sortOrder  String  -SO               10  The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
=========  ======  ========  ==========  ==============================================================================================================================================================

Optional inputs
***************

====================  =======================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
name                  type                     prefix                     position  documentation
====================  =======================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
outputFilename        Optional<Filename>       -O                               10  The sorted SAM/BAM/CRAM output file.
argumentsFile         Optional<Array<File>>    --arguments_file                 10  read one or more arguments files and add them to the command line
compressionLevel      Optional<Integer>        --COMPRESSION_LEVEL              11  Compression level for all compressed files created (e.g. BAM and GELI).
createIndex           Optional<Boolean>        --CREATE_INDEX                   11  Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File         Optional<Boolean>        --CREATE_MD5_FILE                11  Whether to create an MD5 digest for any BAM or FASTQ files created.
maxRecordsInRam       Optional<Integer>        --MAX_RECORDS_IN_RAM             11  When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
quiet                 Optional<Boolean>        --QUIET                          11  Whether to suppress job-summary info on System.err.
reference             Optional<FastaWithDict>  --reference                      11  Reference sequence file.
tmpDir                Optional<String>         --TMP_DIR                        11  Undocumented option
useJdkDeflater        Optional<Boolean>        --use_jdk_deflater               11  Whether to use the JdkDeflater (as opposed to IntelDeflater)
useJdkInflater        Optional<Boolean>        --use_jdk_inflater               11  Whether to use the JdkInflater (as opposed to IntelInflater)
validationStringency  Optional<String>         --VALIDATION_STRINGENCY          11  Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
verbosity             Optional<String>         --verbosity                      11  The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
====================  =======================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: Michael Franklin


*GATK4: SortSAM was last updated on 2019-01-24*.
*This page was automatically generated on 2019-08-02*.
