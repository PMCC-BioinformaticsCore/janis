
.. include:: gatk4_somaticvariantcaller_v0.1.0

GATK4 Somatic Variant Caller
=========================================================

Description
-------------

Tool identifier: ``GATK4_SomaticVariantCaller``

Tool path: ``janis_bioinformatics.tools.variantcallers.gatk.gatksomatic_variants_4_0_12 import GatkSomaticVariantCaller_4_0_12``

Version: v0.1.0





Documentation
-------------

URL
******
*No URL to the documentation was provided*

Tool documentation
******************
This is a VariantCaller based on the GATK Best Practice pipelines. It uses the GATK4 toolkit, specifically 4.0.12.0.

        It has the following steps:

        1. Base Recalibrator x 2
        3. Mutect2
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

===========  ====================  ========  ==========  ===============
name         type                  prefix    position    documentation
===========  ====================  ========  ==========  ===============
normalBam    BamPair
tumorBam     BamPair
normalName   String
tumorName    String
reference    FastaWithDict
snps_dbsnp   CompressedIndexedVCF
snps_1000gp  CompressedIndexedVCF
knownIndels  CompressedIndexedVCF
millsIndels  CompressedIndexedVCF
===========  ====================  ========  ==========  ===============

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


*GATK4 Somatic Variant Caller was last updated on 2019-09-13*.
*This page was automatically generated on 2019-09-26*.
