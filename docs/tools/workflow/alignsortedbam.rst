
Align sorted bam
=================================
*workflow* (alignsortedbam)

Documentation
-------------


URL
******
*No URL to the documentation was provided*: `contribute one <https://github.com/illusional>`_

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
======================  =========  ========  ==========  ===============
name                    type       prefix    position    documentation
======================  =========  ========  ==========  ===============
read_group_header_line  String
fastq                   Fastq
reference               Fasta
tmpdir                  Directory
======================  =========  ========  ==========  ===============

*Align sorted bam was last updated on **Unknown***.

*This page was automatically generated on 2019-01-25*.
