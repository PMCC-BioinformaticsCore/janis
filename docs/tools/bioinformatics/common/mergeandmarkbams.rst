:orphan:

Merge and Mark Duplicates
============================================

*0 contributors · 2 versions*

:ID: ``mergeAndMarkBams``
:Python: ``janis_bioinformatics.tools.common.mergeandmark.mergeandmark_4_1_3 import MergeAndMarkBams_4_1_3``
:Versions: 4.0.12, 4.1.3
:Authors: 
:Citations: 
:Created: None
:Updated: None
:Required inputs:
   - ``bams: Array<IndexedBam>``
:Outputs: 
   - ``out: IndexedBam``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Embedded Tools
***************

======================  ===============================
GATK4: Merge SAM Files  ``Gatk4MergeSamFiles/4.1.3.0``
GATK4: Mark Duplicates  ``Gatk4MarkDuplicates/4.1.3.0``
======================  ===============================

------

Additional configuration (inputs)
---------------------------------

==================================  =================  ================================================================================================================================================================================================================================================================================================================================================================================================
name                                type               documentation
==================================  =================  ================================================================================================================================================================================================================================================================================================================================================================================================
bams                                Array<IndexedBam>
createIndex                         Optional<Boolean>
maxRecordsInRam                     Optional<Integer>
mergeSamFiles_useThreading          Optional<Boolean>  Option to create a background thread to encode, compress and write to disk the output file. The threaded version uses about 20% more CPU and decreases runtime by ~20% when writing out a compressed BAM file.
mergeSamFiles_validationStringency  Optional<String>   Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
==================================  =================  ================================================================================================================================================================================================================================================================================================================================================================================================

.
