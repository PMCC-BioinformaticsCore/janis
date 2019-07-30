
.. include:: gatk4gathervcfs_4.0.12.0

GATK4: Gather VCFs
====================================

Description
-------------

Tool identifier: ``Gatk4GatherVcfs``

Tool path: ``janis_bioinformatics.tools.gatk4.gathervcfs.gathervcfs_4_0 import Gatk4GatherVcfs_4_0``

Version: 4.0.12.0

Container: ``broadinstitute/gatk:4.0.12.0``



Documentation
-------------

URL
******
`https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.12.0/picard_vcf_GatherVcfs.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.12.0/picard_vcf_GatherVcfs.php>`_

Tool documentation
******************
GatherVcfs (Picard)
            
Gathers multiple VCF files from a scatter operation into a single VCF file. 
Input files must be supplied in genomic order and must not have events at overlapping positions.

Outputs
-------
======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======  ==========  ========  ==========  =====================================
name    type        prefix    position    documentation
======  ==========  ========  ==========  =====================================
vcfs    Array<VCF>  --INPUT               [default: []] (-I) Input VCF file(s).
======  ==========  ========  ==========  =====================================

Optional inputs
***************

====================  =====================  =======================  ==========  ======================================================================================================================================================================================================================================================================
name                  type                   prefix                   position    documentation
====================  =====================  =======================  ==========  ======================================================================================================================================================================================================================================================================
outputFilename        Optional<Filename>     --OUTPUT                             [default: null] (-O) Output VCF file.
argumentsFile         Optional<Array<File>>  --arguments_file                     [default: []] read one or more arguments files and add them to the command line
compressionLevel      Optional<Integer>      --COMPRESSION_LEVEL                  [default: 5] Compression level for all compressed files created (e.g. BAM and VCF).
createIndex           Optional<Boolean>      --CREATE_INDEX                       [default: TRUE] Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File         Optional<Boolean>      --CREATE_MD5_FILE                    [default: FALSE] Whether to create an MD5 digest for any BAM or FASTQ files created.
ga4ghClientSecrets    Optional<File>         --GA4GH_CLIENT_SECRETS               [default: client_secrets.json] Google Genomics API client_secrets.json file path.
maxRecordsInRam       Optional<Integer>      --MAX_RECORDS_IN_RAM                 [default: 500000] When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.
quiet                 Optional<Boolean>      --QUIET                              [default: FALSE] Whether to suppress job-summary info on System.err.
referenceSequence     Optional<File>         --REFERENCE_SEQUENCE                 [default: null] Reference sequence file.
tmpDir                Optional<String>       --TMP_DIR                            [default: []] One or more directories with space available to be used by this program for temporary storage of working files
useJdkDeflater        Optional<Boolean>      --USE_JDK_DEFLATER                   [default: FALSE] (-use_jdk_deflater) Use the JDK Deflater instead of the Intel Deflater for writing compressed output
useJdkInflater        Optional<Boolean>      --USE_JDK_INFLATER                   [default: FALSE] (-use_jdk_inflater) Use the JDK Inflater instead of the Intel Inflater for reading compressed input
validationStringency  Optional<String>       --VALIDATION_STRINGENCY              [default: STRICT] Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
verbosity             Optional<Boolean>      --VERBOSITY                          [default: INFO] Control verbosity of logging.
====================  =====================  =======================  ==========  ======================================================================================================================================================================================================================================================================


Metadata
********

Author: Michael Franklin


*GATK4: Gather VCFs was last updated on 2019-05-01*.
*This page was automatically generated on 2019-07-30*.
