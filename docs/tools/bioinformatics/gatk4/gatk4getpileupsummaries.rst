
.. include:: gatk4getpileupsummaries_4.1.3.0
.. include:: gatk4getpileupsummaries_4.1.2.0

GATK4: GetPileupSummaries
===================================================

Description
-------------

Tool identifier: ``GATK4GetPileupSummaries``

Tool path: ``janis_bioinformatics.tools.gatk4.getpileupsummaries.versions import Gatk4GetPileUpSummaries_4_1_3``

Version: 4.1.3.0

Container: ``broadinstitute/gatk:4.1.3.0``

Versions
*********

- 4.1.3.0 (current)
- `4.1.2.0 <gatk4getpileupsummaries_4.1.2.0.html>`_

Documentation
-------------

URL
******
`https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/org_broadinstitute_hellbender_tools_walkers_contamination_GetPileupSummaries.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/org_broadinstitute_hellbender_tools_walkers_contamination_GetPileupSummaries.php>`_

Tool documentation
******************
Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results can be used with CalculateContamination.
The tool requires a common germline variant sites VCF, e.g. the gnomAD resource, with population allele frequencies (AF) in the INFO field. This resource must contain only biallelic SNPs and can be an eight-column sites-only VCF. The tool ignores the filter status of the sites. See the GATK Resource Bundle for an example human file.

Outputs
-------
======  ========  ================================
name    type      documentation
======  ========  ================================
out     TextFile  Table containing the pileup info
======  ========  ================================

Inputs
------
Find the inputs below

Required inputs
***************

======  ====================  ========  ==========  =======================================
name    type                  prefix      position  documentation
======  ====================  ========  ==========  =======================================
bam     Array<BamPair>        -I                 0  The SAM/BAM/CRAM file containing reads.
sites   CompressedIndexedVCF  -V                    sites of common biallelic variants
======  ====================  ========  ==========  =======================================

Optional inputs
***************

==============  ==============================  ===========  ==========  =============================================================
name            type                            prefix         position  documentation
==============  ==============================  ===========  ==========  =============================================================
intervals       Optional<CompressedIndexedVCF>  --intervals              -L (BASE) One or more genomic intervals over which to operate
pileupTableOut  Optional<Filename>              -O                    1
==============  ==============================  ===========  ==========  =============================================================


Metadata
********

Author: Hollizeck Sebastian


*GATK4: GetPileupSummaries was last updated on 2019-09-09*.
*This page was automatically generated on 2019-09-10*.
