:orphan:

GATK4: Merge SAM or BAM with unmapped BAM file
=======================================================================

*2 contributors · 4 versions*

Merges SAM/BAM file with an unmapped BAM file

Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.mergebamalignment.versions import Gatk4MergeBamAlignment_4_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4mergebamalignment_step",
           Gatk4MergeBamAlignment_4_1_2(
               ubam=None,
               bam=None,
           )
       )
       wf.output("out", source=gatk4mergebamalignment_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4MergeBamAlignment:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4MergeBamAlignment > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam:
       - bam_0.sam
       - bam_1.sam
       ubam: ubam.bam




5. Run Gatk4MergeBamAlignment with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4MergeBamAlignment





Information
------------


:ID: ``Gatk4MergeBamAlignment``
:URL: `https://gatk.broadinstitute.org/hc/en-us/articles/360037225832-MergeBamAlignment-Picard- <https://gatk.broadinstitute.org/hc/en-us/articles/360037225832-MergeBamAlignment-Picard->`_
:Versions: 4.1.4.1, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.2.0
:Authors: Michael Franklin (@illusional), Matthias De Smet(@matthdsm)
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2020-02-26



Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     BAM
======  ======  ===============



Additional configuration (inputs)
---------------------------------

=============================  ==========================  ==================================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
name                           type                        prefix                                position  documentation
=============================  ==========================  ==================================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
ubam                           BAM                         --UNMAPPED_BAM                              10  Original SAM or BAM file of unmapped reads, which must be in queryname order.
bam                            Array<SAM>                  --ALIGNED_BAM                               10  SAM or BAM file(s) with alignment data.
reference                      Optional<FastaWithIndexes>  --REFERENCE_SEQUENCE                        10  Reference sequence file.
outputFilename                 Optional<Filename>          --OUTPUT                                    10  Merged SAM or BAM file to write to.
addMateCigar                   Optional<Boolean>           --ADD_MATE_CIGAR                            11  Adds the mate CIGAR tag (MC)
alignedReadsOnly               Optional<Boolean>           --ALIGNED_READS_ONLY                        11  Whether to output only aligned reads.
alignerProperPairFlags         Optional<Boolean>           --ALIGNER_PROPER_PAIR_FLAGS                 11  Use the aligner's idea of what a proper pair is rather than computing in this program.
argumentsFile                  Optional<Array<File>>       --arguments_file                            11  read one or more arguments files and add them to the command line
attributesToRemove             Optional<Array<String>>     --ATTRIBUTES_TO_REMOVE                      11  Attributes from the alignment record that should be removed when merging.
attributesToRetain             Optional<Array<String>>     --ATTRIBUTES_TO_RETAIN                      11  Reserved alignment attributes (tags starting with X, Y, or Z) that should be brought over from the alignment data when merging.
attributesToReverse            Optional<Array<String>>     --ATTRIBUTES_TO_REVERSE                     11  Attributes on negative strand reads that need to be reversed.
attributesToReverseComplement  Optional<Array<String>>     --ATTRIBUTES_TO_REVERSE_COMPLEMENT          11  Attributes on negative strand reads that need to be reverse complemented.
clipAdapter                    Optional<Boolean>           --CLIP_ADAPTERS                             11  Whether to clip adapters where identified.
clipOverlappingReads           Optional<Boolean>           --CLIP_OVERLAPPING_READS                    11  For paired reads, soft clip the 3' end of each read if necessary so that it does not extend past the 5' end of its mate.
expectedOrientations           Optional<Array<String>>     --EXPECTED_ORIENTATIONS                     11  The expected orientation of proper read pairs.
includeSecondaryAlginments     Optional<Boolean>           --INCLUDE_SECONDARY_ALIGNMENTS              11  If false, do not write secondary alignments to output.
isBisulfiteSequencing          Optional<Boolean>           --IS_BISULFITE_SEQUENCE                     11  Whether the lane is bisulfite sequence (used when calculating the NM tag).
matchingDictionaryTags         Optional<Array<String>>     --MATCHING_DICTIONARY_TAGS                  11  List of Sequence Records tags that must be equal (if present) in the reference dictionary and in the aligned file.
maxInsertionsOrDeletions       Optional<Integer>           --MAX_INSERTIONS_OR_DELETIONS               11  The maximum number of insertions or deletions permitted for an alignment to be included.
minUnclippedBases              Optional<Integer>           --MIN_UNCLIPPED_BASES                       11  If UNMAP_CONTAMINANT_READS is set, require this many unclipped bases or else the read will be marked as contaminant.
primaryAlignmentStrategy       Optional<Integer>           --PRIMARY_ALIGNMENT_STRATEGY                11  Strategy for selecting primary alignment when the aligner has provided more than one alignment for a pair or fragment, and none are marked as primary, more than one is marked as primary, or the primary alignment is filtered out for some reason.
programGroupCommandLine        Optional<String>            --PROGRAM_GROUP_COMMAND_LINE                11  The command line of the program group.
programGroupName               Optional<String>            --PROGRAM_GROUP_NAME                        11  The name of the program group.
programGroupVersion            Optional<String>            --PROGRAM_GROUP_VERSION                     11  The version of the program group.
programRecordId                Optional<String>            --PROGRAM_RECORD_ID                         11  The program group ID of the aligner.
sortOrder                      Optional<String>            -SO                                         10  The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
unmapContaminantReads          Optional<Boolean>           --UNMAP_CONTAMINANT_READS                   11  Detect reads originating from foreign organisms (e.g. bacterial DNA in a non-bacterial sample),and unmap + label those reads accordingly.
unmappedReadStrategy           Optional<String>            --UNMAPPED_READ_STRATEGY                    11  How to deal with alignment information in reads that are being unmapped (e.g. due to cross-species contamination.) Currently ignored unless UNMAP_CONTAMINANT_READS = true.
addPgTagToReads                Optional<Boolean>           --ADD_PG_TAG_TO_READS                       11  Add PG tag to each read in a SAM or BAM
compressionLevel               Optional<Integer>           --COMPRESSION_LEVEL                         11  Compression level for all compressed files created (e.g. BAM and GELI).
createIndex                    Optional<Boolean>           --CREATE_INDEX                              11  Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File                  Optional<Boolean>           --CREATE_MD5_FILE                           11  Whether to create an MD5 digest for any BAM or FASTQ files created.
maxRecordsInRam                Optional<Integer>           --MAX_RECORDS_IN_RAM                        11  When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
quiet                          Optional<Boolean>           --QUIET                                     11  Whether to suppress job-summary info on System.err.
tmpDir                         Optional<String>            --TMP_DIR                                   11  Undocumented option
useJdkDeflater                 Optional<Boolean>           --use_jdk_deflater                          11  Whether to use the JdkDeflater (as opposed to IntelDeflater)
useJdkInflater                 Optional<Boolean>           --use_jdk_inflater                          11  Whether to use the JdkInflater (as opposed to IntelInflater)
validationStringency           Optional<String>            --VALIDATION_STRINGENCY                     11  Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
verbosity                      Optional<String>            --verbosity                                 11  The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
=============================  ==========================  ==================================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
