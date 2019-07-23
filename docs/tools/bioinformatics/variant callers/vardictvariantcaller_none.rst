:orphan:


Vardict Variant Caller
=============================================
Tool identifier: ``vardictVariantCaller``
Tool path: ``janis_bioinformatics.tools.variantcallers.vardictgermline_variants import VardictGermlineVariantCaller``

Version: None




Documentation
-------------

URL
******
*No URL to the documentation was provided*

Description
*********
*No documentation was provided: `contribute one <https://github.com/illusional>`_*

Outputs
-------
===============  ======  ===============
name             type    documentation
===============  ======  ===============
vardictVariants  VCF
out              VCF
===============  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

==================  =============  ========  ==========  ===============
name                type           prefix    position    documentation
==================  =============  ========  ==========  ===============
intervals           bed
bam                 BamPair
reference           FastaWithDict
sampleName          String
allelFreqThreshold  Float
headerLines         File
==================  =============  ========  ==========  ===============

Optional inputs
***************

====================  =================  ========  ==========  ===============
name                  type               prefix    position    documentation
====================  =================  ========  ==========  ===============
chromNamesAreNumbers  Optional<Boolean>
vcfFormat             Optional<Boolean>
chromColumn           Optional<Integer>
regStartCol           Optional<Integer>
geneEndCol            Optional<Integer>
====================  =================  ========  ==========  ===============


Metadata
********

Author: **Unknown**


*Vardict Variant Caller was last updated on **Unknown***.
*This page was automatically generated on 2019-07-23*.
