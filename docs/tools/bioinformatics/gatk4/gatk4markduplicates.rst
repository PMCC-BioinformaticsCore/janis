
.. include:: gatk4markduplicates_4.0.12.0

GATK4: Mark Duplicates
============================================

Description
-------------

Tool identifier: ``Gatk4MarkDuplicates``

Tool path: ``janis_bioinformatics.tools.gatk4.markduplicates.latest import Gatk4MarkDuplicatesLatest``

Version: 4.0.12.0

Docker: ``broadinstitute/gatk:4.0.12.0``



Documentation
-------------

URL
******
`https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_markduplicates_MarkDuplicates.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_markduplicates_MarkDuplicates.php>`_

Description
*********
MarkDuplicates (Picard): Identifies duplicate reads.

This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are 
defined as originating from a single fragment of DNA. Duplicates can arise during sample 
preparation e.g. library construction using PCR. See also EstimateLibraryComplexity for 
additional notes on PCR duplication artifacts. Duplicate reads can also result from a single 
amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the 
sequencing instrument. These duplication artifacts are referred to as optical duplicates.

The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads 
and read-pairs in a SAM/BAM file. An BARCODE_TAG option is available to facilitate duplicate
marking using molecular barcodes. After duplicate reads are collected, the tool differentiates 
the primary and duplicate reads using an algorithm that ranks reads by the sums of their 
base-quality scores (default method).

The tool's main output is a new SAM or BAM file, in which duplicates have been identified 
in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, 
which corresponds to a decimal value of 1024. If you are not familiar with this type of annotation, 
please see the following blog post for additional information.

Although the bitwise flag annotation indicates whether a read was marked as a duplicate, 
it does not identify the type of duplicate. To do this, a new tag called the duplicate type (DT) 
tag was recently added as an optional output in the 'optional field' section of a SAM/BAM file. 
Invoking the TAGGING_POLICY option, you can instruct the program to mark all the duplicates (All), 
only the optical duplicates (OpticalOnly), or no duplicates (DontTag). The records within the 
output of a SAM/BAM file will have values for the 'DT' tag (depending on the invoked TAGGING_POLICY), 
as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ). 
This tool uses the READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the 
primary methods to identify and differentiate duplicate types. Set READ_NAME_REGEX to null to 
skip optical duplicate detection, e.g. for RNA-seq or other data where duplicate sets are 
extremely large and estimating library complexity is not an aim. Note that without optical 
duplicate counts, library size estimation will be inaccurate.

MarkDuplicates also produces a metrics file indicating the numbers 
of duplicates for both single- and paired-end reads.

The program can take either coordinate-sorted or query-sorted inputs, however the behavior 
is slightly different. When the input is coordinate-sorted, unmapped mates of mapped records 
and supplementary/secondary alignments are not marked as duplicates. However, when the input 
is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary 
reads are not excluded from the duplication test and can be marked as duplicate reads.

If desired, duplicates can be removed using the REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES options.

Outputs
-------
=======  =======  ===============
name     type     documentation
=======  =======  ===============
out      BamPair
metrics  tsv
=======  =======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======  =======  ========  ==========  =========================================================================
name    type     prefix      position  documentation
======  =======  ========  ==========  =========================================================================
bam     BamPair  -I                10  One or more input SAM or BAM files to analyze. Must be coordinate sorted.
======  =======  ========  ==========  =========================================================================

Optional inputs
***************

====================  =======================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
name                  type                     prefix                     position  documentation
====================  =======================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
outputFilename        Optional<Filename>       -O                               10  File to write duplication metrics to
metricsFilename       Optional<Filename>       -M                               10  The output file to write marked records to.
argumentsFile         Optional<Array<File>>    --arguments_file                 10  read one or more arguments files and add them to the command line
assumeSortOrder       Optional<String>         -ASO                                 If not null, assume that the input file has this order even if the header says otherwise. Exclusion: This argument cannot be used at the same time as ASSUME_SORTED. The --ASSUME_SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
barcodeTag            Optional<String>         --BARCODE_TAG                        Barcode SAM tag (ex. BC for 10X Genomics)
comment               Optional<Array<String>>  -CO                                  Comment(s) to include in the output file's header.
compressionLevel      Optional<Integer>        --COMPRESSION_LEVEL              11  Compression level for all compressed files created (e.g. BAM and GELI).
createIndex           Optional<Boolean>        --CREATE_INDEX                   11  Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File         Optional<Boolean>        --CREATE_MD5_FILE                11  Whether to create an MD5 digest for any BAM or FASTQ files created.
maxRecordsInRam       Optional<Integer>        --MAX_RECORDS_IN_RAM             11  When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
quiet                 Optional<Boolean>        --QUIET                          11  Whether to suppress job-summary info on System.err.
tmpDir                Optional<String>         --TMP_DIR                        11  Undocumented option
useJdkDeflater        Optional<Boolean>        --use_jdk_deflater               11  Whether to use the JdkDeflater (as opposed to IntelDeflater)
useJdkInflater        Optional<Boolean>        --use_jdk_inflater               11  Whether to use the JdkInflater (as opposed to IntelInflater)
validationStringency  Optional<String>         --VALIDATION_STRINGENCY          11  Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
verbosity             Optional<String>         --verbosity                      11  The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
====================  =======================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: Michael Franklin


*GATK4: Mark Duplicates was last updated on 2019-01-24*.
*This page was automatically generated on 2019-07-24*.
