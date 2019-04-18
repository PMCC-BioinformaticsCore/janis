
Align sorted BAM
=================================
Tool identifier: ``alignsortedbam``

Tool path: ``from janis_bioinformatics.tools import AlignSortedBam``

Documentation
-------------


URL
******
*No URL to the documentation was provided*

Docstring
*********
Align sorted bam with this subworkflow consisting of BWA Mem + SamTools + Gatk4SortSam

Outputs
-------
============  =======  ===============
name          type     documentation
============  =======  ===============
out_bwa       SAM
out_samtools  BAM
out           BamPair
============  =======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

===================  =============  ========  ==========  ===============
name                 type           prefix    position    documentation
===================  =============  ========  ==========  ===============
fastq                Fastq
readGroupHeaderLine  String
reference            FastaWithDict
===================  =============  ========  ==========  ===============

Optional inputs
***************

====================  =================  ========  ==========  ===============
name                  type               prefix    position    documentation
====================  =================  ========  ==========  ===============
adapter               Optional<String>
adapter_g             Optional<String>
removeMiddle5Adapter  Optional<String>
removeMiddle3Adapter  Optional<String>
qualityCutoff         Optional<Integer>
minReadLength         Optional<Integer>
sortOrder             Optional<String>
createIndex           Optional<Boolean>
validationStringency  Optional<String>
maxRecordsInRam       Optional<Integer>
====================  =================  ========  ==========  ===============


Metadata
********

Author: Michael Franklin


*Align sorted BAM was last updated on **Unknown***.
*This page was automatically generated on 2019-04-18*.
