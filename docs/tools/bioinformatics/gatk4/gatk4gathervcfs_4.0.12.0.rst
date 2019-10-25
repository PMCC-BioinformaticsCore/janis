:orphan:

GATK4: Gather VCFs
====================================

1 contributor · 4 versions

:ID: ``Gatk4GatherVcfs``
:Python: ``janis_bioinformatics.tools.gatk4.gathervcfs.versions import Gatk4GatherVcfs_4_0``
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.0.12.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-05-01
:Updated: 2019-05-01
:Required inputs:
   - ``vcfs: Array<VCF>``
:Outputs: 
   - ``out: VCF``

Documentation
-------------

URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.12.0/picard_vcf_GatherVcfs.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.12.0/picard_vcf_GatherVcfs.php>`_

GatherVcfs (Picard)
            
Gathers multiple VCF files from a scatter operation into a single VCF file. 
Input files must be supplied in genomic order and must not have events at overlapping positions.

------

Additional configuration (inputs)
---------------------------------

====================  =====================  ======================================================================================================================================================================================================================================================================
name                  type                   documentation
====================  =====================  ======================================================================================================================================================================================================================================================================
vcfs                  Array<VCF>             [default: []] (-I) Input VCF file(s).
outputFilename        Optional<Filename>     [default: null] (-O) Output VCF file.
argumentsFile         Optional<Array<File>>  [default: []] read one or more arguments files and add them to the command line
compressionLevel      Optional<Integer>      [default: 5] Compression level for all compressed files created (e.g. BAM and VCF).
createIndex           Optional<Boolean>      [default: TRUE] Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File         Optional<Boolean>      [default: FALSE] Whether to create an MD5 digest for any BAM or FASTQ files created.
ga4ghClientSecrets    Optional<File>         [default: client_secrets.json] Google Genomics API client_secrets.json file path.
maxRecordsInRam       Optional<Integer>      [default: 500000] When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.
quiet                 Optional<Boolean>      [default: FALSE] Whether to suppress job-summary info on System.err.
referenceSequence     Optional<File>         [default: null] Reference sequence file.
tmpDir                Optional<String>       [default: []] One or more directories with space available to be used by this program for temporary storage of working files
useJdkDeflater        Optional<Boolean>      [default: FALSE] (-use_jdk_deflater) Use the JDK Deflater instead of the Intel Deflater for writing compressed output
useJdkInflater        Optional<Boolean>      [default: FALSE] (-use_jdk_inflater) Use the JDK Inflater instead of the Intel Inflater for reading compressed input
validationStringency  Optional<String>       [default: STRICT] Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
verbosity             Optional<Boolean>      [default: INFO] Control verbosity of logging.
====================  =====================  ======================================================================================================================================================================================================================================================================

