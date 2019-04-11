
Align sorted BAM
=================================
Tool identifier: ``alignsortedbam``

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

======================  =============  ========  ==========  ===============
name                    type           prefix    position    documentation
======================  =============  ========  ==========  ===============
fastq                   Fastq
read_group_header_line  String
reference               FastaWithDict
======================  =============  ========  ==========  ===============

Optional inputs
***************

======  ======  ========  ==========  ===============
name    type    prefix    position    documentation
======  ======  ========  ==========  ===============
======  ======  ========  ==========  ===============


Metadata
********

Author: Michael Franklin


*Align sorted BAM was last updated on **Unknown***.
*This page was automatically generated on 2019-04-11*.
