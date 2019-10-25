:orphan:

GATK4: SortSAM
=============================

1 contributor · 4 versions

:ID: ``gatk4sortsam``
:Python: ``janis_bioinformatics.tools.gatk4.sortsam.versions import Gatk4SortSam_4_1_3``
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.3.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24
:Required inputs:
   - ``bam: BAM``

   - ``sortOrder: String``
:Outputs: 
   - ``out: BamPair``

Documentation
-------------

URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.3/org_broadinstitute_hellbender_tools_picard_sam_SortSam.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.3/org_broadinstitute_hellbender_tools_picard_sam_SortSam.php>`_

Sorts a SAM/BAM/CRAM file.

------

Additional configuration (inputs)
---------------------------------

====================  =======================  ================================================================================================================================================================================================================================================================================================================================================================================================
name                  type                     documentation
====================  =======================  ================================================================================================================================================================================================================================================================================================================================================================================================
bam                   BAM                      The SAM/BAM/CRAM file to sort.
sortOrder             String                   The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
outputFilename        Optional<Filename>       The sorted SAM/BAM/CRAM output file.
argumentsFile         Optional<Array<File>>    read one or more arguments files and add them to the command line
compressionLevel      Optional<Integer>        Compression level for all compressed files created (e.g. BAM and GELI).
createIndex           Optional<Boolean>        Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File         Optional<Boolean>        Whether to create an MD5 digest for any BAM or FASTQ files created.
maxRecordsInRam       Optional<Integer>        When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
quiet                 Optional<Boolean>        Whether to suppress job-summary info on System.err.
reference             Optional<FastaWithDict>  Reference sequence file.
tmpDir                Optional<String>         Undocumented option
useJdkDeflater        Optional<Boolean>        Whether to use the JdkDeflater (as opposed to IntelDeflater)
useJdkInflater        Optional<Boolean>        Whether to use the JdkInflater (as opposed to IntelInflater)
validationStringency  Optional<String>         Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
verbosity             Optional<String>         The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
====================  =======================  ================================================================================================================================================================================================================================================================================================================================================================================================

