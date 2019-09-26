:orphan:


Strelka Germline Variant Caller
==============================================================

Description
-------------

Tool identifier: ``strelkaGermlineVariantCaller``

Tool path: ``janis_bioinformatics.tools.variantcallers.illuminagermline_strelka import IlluminaGermlineVariantCaller``

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
========  ====================  ===============
name      type                  documentation
========  ====================  ===============
diploid   CompressedIndexedVCF
variants  CompressedIndexedVCF
out       VCF
========  ====================  ===============

Inputs
------
Find the inputs below

Required inputs
***************

=========  =============  ========  ==========  ===============
name       type           prefix    position    documentation
=========  =============  ========  ==========  ===============
bam        BamPair
reference  FastaWithDict
=========  =============  ========  ==========  ===============

Optional inputs
***************

====================  =======================  ========  ==========  ===============
name                  type                     prefix    position    documentation
====================  =======================  ========  ==========  ===============
intervals             Optional<BedTABIX>
bcfview_applyFilters  Optional<Array<String>>
====================  =======================  ========  ==========  ===============


Metadata
********

Author: **Unknown**


*Strelka Germline Variant Caller was last updated on **Unknown***.
*This page was automatically generated on 2019-09-26*.
