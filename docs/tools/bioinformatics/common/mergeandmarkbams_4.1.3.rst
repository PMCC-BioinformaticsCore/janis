:orphan:

Merge and Mark Duplicates
============================================

0 contributors · 2 versions

:ID: ``mergeAndMarkBams``
:Python: ``janis_bioinformatics.tools.common.mergeandmark.mergeandmark_4_1_3 import MergeAndMarkBams_4_1_3``
:Versions: 4.0.12, 4.1.3
:Authors: 
:Citations: 
:Created: None
:Updated: None
:Required inputs:
   - ``bams: Array<BamPair>``
:Outputs: 
   - ``out: BamPair``

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

==================================  =================  ===============
name                                type               documentation
==================================  =================  ===============
bams                                Array<BamPair>
createIndex                         Optional<Boolean>
maxRecordsInRam                     Optional<Integer>
mergeSamFiles_useThreading          Optional<Boolean>
mergeSamFiles_validationStringency  Optional<String>
==================================  =================  ===============

.
