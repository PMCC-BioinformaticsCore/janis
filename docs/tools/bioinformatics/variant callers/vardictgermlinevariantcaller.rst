
.. include:: vardictgermlinevariantcaller_v0.1.0

Vardict Germline Variant Caller
==============================================================

Description
-------------

Tool identifier: ``vardictGermlineVariantCaller``

Tool path: ``janis_bioinformatics.tools.variantcallers.vardictgermline_variants import VardictGermlineVariantCaller``

Version: v0.1.0





Documentation
-------------

URL
******
*No URL to the documentation was provided*

Tool documentation
******************
*No documentation was provided: `contribute one <https://github.com/illusional>`_*

Outputs
-------
===============  =============  ===============
name             type           documentation
===============  =============  ===============
vardictVariants  CompressedVCF
out              VCF
===============  =============  ===============

Inputs
------
Find the inputs below

Required inputs
***************

===========  =============  ========  ==========  ===============
name         type           prefix    position    documentation
===========  =============  ========  ==========  ===============
bam          BamPair
intervals    bed
sampleName   String
headerLines  File
reference    FastaWithDict
===========  =============  ========  ==========  ===============

Optional inputs
***************

============================  =================  ========  ==========  ===============
name                          type               prefix    position    documentation
============================  =================  ========  ==========  ===============
alleleFreqThreshold           Optional<Float>
vardict_chromNamesAreNumbers  Optional<Boolean>
vardict_vcfFormat             Optional<Boolean>
vardict_chromColumn           Optional<Integer>
vardict_regStartCol           Optional<Integer>
vardict_geneEndCol            Optional<Integer>
============================  =================  ========  ==========  ===============


Metadata
********

Author: **Unknown**


*Vardict Germline Variant Caller was last updated on **Unknown***.
*This page was automatically generated on 2019-09-26*.
