
.. include:: alignsortedbam_1.0.0

Align sorted BAM
=================================

Description
-------------

Tool identifier: ``alignsortedbam``

Tool path: ``janis_bioinformatics.tools.common.alignsortedbam import AlignSortedBam``

Version: 1.0.0





Documentation
-------------

URL
******
*No URL to the documentation was provided*

Tool documentation
******************
Align sorted bam with this subworkflow consisting of BWA Mem + SamTools + Gatk4SortSam

Outputs
-------
=======  =======  ===============
name     type     documentation
=======  =======  ===============
out_bwa  BAM
out      BamPair
=======  =======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

==========  =============  ========  ==========  ===============
name        type           prefix    position    documentation
==========  =============  ========  ==========  ===============
fastq       Fastq
sampleName  String
reference   FastaWithDict
==========  =============  ========  ==========  ===============

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
sortSamTmpDir         Optional<String>
====================  =================  ========  ==========  ===============


Metadata
********

Author: Michael Franklin


*Align sorted BAM was last updated on **Unknown***.
*This page was automatically generated on 2019-07-26*.
