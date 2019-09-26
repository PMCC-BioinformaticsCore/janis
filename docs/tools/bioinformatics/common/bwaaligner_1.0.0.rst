:orphan:


Align and sort reads
=================================

Description
-------------

Tool identifier: ``BwaAligner``

Tool path: ``janis_bioinformatics.tools.common.bwaaligner import BwaAligner``

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
======  =======  ===============
name    type     documentation
======  =======  ===============
out     BamPair
======  =======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

=========  =============  ========  ==========  ===============
name       type           prefix    position    documentation
=========  =============  ========  ==========  ===============
name       String
reference  FastaWithDict
fastq      FastqGzPair
=========  =============  ========  ==========  ===============

Optional inputs
***************

=============================  =================  ========  ==========  ===============
name                           type               prefix    position    documentation
=============================  =================  ========  ==========  ===============
cutadapt_adapter               Optional<String>
cutadapt_adapter_g             Optional<String>
cutadapt_removeMiddle5Adapter  Optional<String>
cutadapt_removeMiddle3Adapter  Optional<String>
cutadapt_qualityCutoff         Optional<Integer>
cutadapt_minReadLength         Optional<Integer>
bwamem_markShorterSplits       Optional<Boolean>
sortsam_sortOrder              Optional<String>
sortsam_createIndex            Optional<Boolean>
sortsam_validationStringency   Optional<String>
sortsam_maxRecordsInRam        Optional<Integer>
sortsam_tmpDir                 Optional<String>
=============================  =================  ========  ==========  ===============


Metadata
********

Author: Michael Franklin


*Align and sort reads was last updated on **Unknown***.
*This page was automatically generated on 2019-09-26*.
