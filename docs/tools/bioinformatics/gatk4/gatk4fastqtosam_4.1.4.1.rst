:orphan:

GATK4: Convert a FASTQ file to an unaligned BAM or SAM file.
==============================================================================

*2 contributors · 4 versions*

:ID: ``Gatk4FastqToSam``
:Python: ``janis_bioinformatics.tools.gatk4.fastqtosam.versions import Gatk4FastqToSam_4_1_4``
:Versions: 4.1.4.1, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.4.1
:Authors: Michael Franklin (@illisional), Matthias De Smet(@matthdsm)
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2020-02-26
:Updated: 2020-02-26
:Required inputs:
   - ``fastqR1: FastqGz``
:Outputs: 
   - ``out: IndexedBam``

Documentation
-------------

URL: `https://gatk.broadinstitute.org/hc/en-us/articles/360037226792-FastqToSam-Picard- <https://gatk.broadinstitute.org/hc/en-us/articles/360037226792-FastqToSam-Picard->`_

Converts a FASTQ file to an unaligned BAM or SAM file.

------

None

Additional configuration (inputs)
---------------------------------

========================  ==========================  ==============================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
name                      type                        prefix                            position  documentation
========================  ==========================  ==============================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
fastqR1                   FastqGz                     --FASTQ                                 10  Input fastq file (optionally gzipped) for single end data, or first read in paired end data.
fastqR2                   Optional<FastqGz>           --FASTQ2                                10  Input fastq file (optionally gzipped) for single end data, or first read in paired end data.
sampleName                Optional<String>            --SAMPLE_NAME                           10  Input fastq file (optionally gzipped) for single end data, or first read in paired end data.
reference                 Optional<FastaWithIndexes>  --REFERENCE_SEQUENCE                    10  Reference sequence file.
outputFilename            Optional<Filename>          --OUTPUT                                10  Merged SAM or BAM file to write to.
allowAndIgnoreEmptyLines  Optional<Boolean>           --ALLOW_AND_IGNORE_EMPTY_LINES          11  Allow (and ignore) empty lines
argumentsFile             Optional<Array<File>>       --arguments_file                        11  read one or more arguments files and add them to the command line
comment                   Optional<Array<String>>     --COMMENT                               11  Comment(s) to include in the merged output file's header.
description               Optional<Array<String>>     --DESCRIPTION                           11  Inserted into the read group header
libraryName               Optional<Array<String>>     --LIBRARY_NAME                          11  The library name to place into the LB attribute in the read group header
maxQ                      Optional<Integer>           --MAX_Q                                 11  Maximum quality allowed in the input fastq. An exception will be thrown if a quality is greater than this value.
minQ                      Optional<Integer>           --MIN_Q                                 11  Minimum quality allowed in the input fastq. An exception will be thrown if a quality is less than this value.
platform                  Optional<String>            --PLATFORM                              11  The platform type (e.g. ILLUMINA, SOLID) to insert into the read group header.
platformModel             Optional<String>            --PLATFORM_MODEL                        11  Platform model to insert into the group header (free-form text providing further details of the platform/technology used).
platformUnit              Optional<String>            --PLATFORM_UNIT                         11  The expected orientation of proper read pairs.
predictedInsertSize       Optional<Integer>           --PREDICTED_INSERT_SIZE                 11  Predicted median insert size, to insert into the read group header.
programGroup              Optional<String>            --PROGRAM_GROUP                         11  Program group to insert into the read group header.
readGroupName             Optional<String>            --READ_GROUP_NAME                       11  Read group name.
runDate                   Optional<String>            --RUN_DATE                              11  Date the run was produced, to insert into the read group header
sequencingCenter          Optional<String>            --SEQUENCING_CENTER                     11  The sequencing center from which the data originated.
sortOrder                 Optional<String>            -SO                                     10  The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
useSequenctialFastqs      Optional<Boolean>           --USE_SEQUENTIAL_FASTQS                 11  Use sequential fastq files with the suffix _###.fastq or _###.fastq.gz.
compressionLevel          Optional<Integer>           --COMPRESSION_LEVEL                     11  Compression level for all compressed files created (e.g. BAM and GELI).
createIndex               Optional<Boolean>           --CREATE_INDEX                          11  Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File             Optional<Boolean>           --CREATE_MD5_FILE                       11  Whether to create an MD5 digest for any BAM or FASTQ files created.
maxRecordsInRam           Optional<Integer>           --MAX_RECORDS_IN_RAM                    11  When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
quiet                     Optional<Boolean>           --QUIET                                 11  Whether to suppress job-summary info on System.err.
tmpDir                    Optional<String>            --TMP_DIR                               11  Undocumented option
useJdkDeflater            Optional<Boolean>           --use_jdk_deflater                      11  Whether to use the JdkDeflater (as opposed to IntelDeflater)
useJdkInflater            Optional<Boolean>           --use_jdk_inflater                      11  Whether to use the JdkInflater (as opposed to IntelInflater)
validationStringency      Optional<String>            --VALIDATION_STRINGENCY                 11  Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
verbosity                 Optional<String>            --verbosity                             11  The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
========================  ==========================  ==============================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================

