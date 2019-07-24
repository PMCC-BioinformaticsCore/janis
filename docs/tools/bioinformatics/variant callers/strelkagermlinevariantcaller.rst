
.. include:: strelkagermlinevariantcaller_none

Strelka Germline Variant Caller
==============================================================

Description
-------------

Tool identifier: ``strelkaGermlineVariantCaller``

Tool path: ``janis_bioinformatics.tools.variantcallers.illuminagermline_strelka import IlluminaGermlineVariantCaller``

Version: None





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
========  ==========  ===============
name      type        documentation
========  ==========  ===============
diploid   vcf-gz-tbi
variants  vcf-gz-tbi
out       VCF
========  ==========  ===============

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

==============  =======================  ========  ==========  ===============
name            type                     prefix    position    documentation
==============  =======================  ========  ==========  ===============
strelkaRegions  Optional<BedTABIX>
filters         Optional<Array<String>>
==============  =======================  ========  ==========  ===============


Metadata
********

Author: **Unknown**


*Strelka Germline Variant Caller was last updated on **Unknown***.
*This page was automatically generated on 2019-07-24*.
