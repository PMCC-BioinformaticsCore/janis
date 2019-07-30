
.. include:: gatk4_germlinevariantcaller_0.1.0%22

GATK4 Germline Variant Caller
===========================================================

Description
-------------

Tool identifier: ``GATK4_GermlineVariantCaller``

Tool path: ``janis_bioinformatics.tools.variantcallers.gatkgermline_variants import GatkGermlineVariantCaller``

Version: 0.1.0"





Documentation
-------------

URL
******
*No URL to the documentation was provided*

Tool documentation
******************
This is a VariantCaller based on the GATK Best Practice pipelines. It uses the GATK4 toolkit, specifically 4.0.12.0.

It has the following steps:

1. BaseRecalibrator
2. ApplyBQSR
3. HaplotypeCaller
4. SplitMultiAllele

Outputs
-------
======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

===========  =============  ========  ==========  ===============
name         type           prefix    position    documentation
===========  =============  ========  ==========  ===============
bam          BamPair
reference    FastaWithDict
snps_dbsnp   vcf-gz-tbi
snps_1000gp  vcf-gz-tbi
knownIndels  vcf-gz-tbi
millsIndels  vcf-gz-tbi
===========  =============  ========  ==========  ===============

Optional inputs
***************

=========  =============  ========  ==========  ===============
name       type           prefix    position    documentation
=========  =============  ========  ==========  ===============
intervals  Optional<bed>
=========  =============  ========  ==========  ===============


Metadata
********

Author: **Unknown**


*GATK4 Germline Variant Caller was last updated on **Unknown***.
*This page was automatically generated on 2019-07-30*.
