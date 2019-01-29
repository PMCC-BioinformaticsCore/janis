
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
===========  =======  ===============
name         type     documentation
===========  =======  ===============
o1_bwa       SAM
o2_samtools  BAM
o3_sortsam   BamPair
===========  =======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======================  =========  ========  ==========  ===============
name                    type       prefix    position    documentation
======================  =========  ========  ==========  ===============
read_group_header_line  String
fastq                   Fastq
reference               Fasta
tmpdir                  Directory
======================  =========  ========  ==========  ===============

Optional inputs
***************

======  ======  ========  ==========  ===============
name    type    prefix    position    documentation
======  ======  ========  ==========  ===============
======  ======  ========  ==========  ===============


*Align sorted BAM was last updated on **Unknown***.
*This page was automatically generated on 2019-01-30*.
