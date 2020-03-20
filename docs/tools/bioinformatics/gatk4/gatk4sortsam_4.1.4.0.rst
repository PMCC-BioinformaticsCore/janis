:orphan:

GATK4: SortSAM
=============================

*1 contributor · 4 versions*

Sorts a SAM/BAM/CRAM file.

Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.sortsam.versions import Gatk4SortSam_4_1_4

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4sortsam_step",
           Gatk4SortSam_4_1_4(
               bam=None,
               sortOrder=None,
           )
       )
       wf.output("out", source=gatk4sortsam_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4SortSam:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4SortSam > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       sortOrder: <value>




5. Run Gatk4SortSam with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4SortSam





Information
------------


:ID: ``Gatk4SortSam``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.3/org_broadinstitute_hellbender_tools_picard_sam_SortSam.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.3/org_broadinstitute_hellbender_tools_picard_sam_SortSam.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.4.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24



Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedBam
======  ==========  ===============



Additional configuration (inputs)
---------------------------------

====================  ==========================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
name                  type                        prefix                     position  documentation
====================  ==========================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
bam                   BAM                         -I                               10  The SAM/BAM/CRAM file to sort.
sortOrder             String                      -SO                              10  The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
outputFilename        Optional<Filename>          -O                               10  The sorted SAM/BAM/CRAM output file.
argumentsFile         Optional<Array<File>>       --arguments_file                 10  read one or more arguments files and add them to the command line
compressionLevel      Optional<Integer>           --COMPRESSION_LEVEL              11  Compression level for all compressed files created (e.g. BAM and GELI).
createIndex           Optional<Boolean>           --CREATE_INDEX                   11  Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File         Optional<Boolean>           --CREATE_MD5_FILE                11  Whether to create an MD5 digest for any BAM or FASTQ files created.
maxRecordsInRam       Optional<Integer>           --MAX_RECORDS_IN_RAM             11  When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
quiet                 Optional<Boolean>           --QUIET                          11  Whether to suppress job-summary info on System.err.
reference             Optional<FastaWithIndexes>  --reference                      11  Reference sequence file.
tmpDir                Optional<String>            --TMP_DIR                        11  Undocumented option
useJdkDeflater        Optional<Boolean>           --use_jdk_deflater               11  Whether to use the JdkDeflater (as opposed to IntelDeflater)
useJdkInflater        Optional<Boolean>           --use_jdk_inflater               11  Whether to use the JdkInflater (as opposed to IntelInflater)
validationStringency  Optional<String>            --VALIDATION_STRINGENCY          11  Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
verbosity             Optional<String>            --verbosity                      11  The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
====================  ==========================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
