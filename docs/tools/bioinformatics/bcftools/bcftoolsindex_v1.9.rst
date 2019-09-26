:orphan:


BCFTools: Index
===============================

Description
-------------

Tool identifier: ``bcftoolsIndex``

Tool path: ``janis_bioinformatics.tools.bcftools.index.versions import BcfToolsIndex_1_9``

Version: v1.9

Container: ``michaelfranklin/bcftools:1.9``



Documentation
-------------

URL
******
`https://samtools.github.io/bcftools/bcftools.html#norm <https://samtools.github.io/bcftools/bcftools.html#norm>`_

Tool documentation
******************
Index bgzip compressed VCF/BCF files for random access.

Outputs
-------
======  ====================  ===============
name    type                  documentation
======  ====================  ===============
out     CompressedIndexedVCF
======  ====================  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======  =============  ========  ==========  ===============
name    type           prefix      position  documentation
======  =============  ========  ==========  ===============
vcf     CompressedVCF                     1
======  =============  ========  ==========  ===============

Optional inputs
***************

========  =================  ===========  ==========  ============================================================
name      type               prefix       position    documentation
========  =================  ===========  ==========  ============================================================
csi       Optional<Boolean>  --csi                    (-c) generate CSI-format index for VCF/BCF files [default]
force     Optional<Boolean>  --force                  (-f) overwrite index if it already exists
minShift  Optional<Integer>  --min-shift              (-m) set minimal interval size for CSI indices to 2^INT [14]
tbi       Optional<Boolean>  --tbi                    (-t) generate TBI-format index for VCF files
threads   Optional<Integer>  --threads                sets the number of threads [0]
nrecords  Optional<Boolean>  --nrecords               (-n) print number of records based on existing index file
stats     Optional<Boolean>  --stats                  (-s) print per contig stats based on existing index file
========  =================  ===========  ==========  ============================================================


Metadata
********

Author: **Unknown**


*BCFTools: Index was last updated on **Unknown***.
*This page was automatically generated on 2019-09-26*.
