
Strelka Variant Caller
=============================================
Tool identifier: ``strelkaVariantCaller``

Tool path: ``from janis_bioinformatics.tools.illumina import StrelkaVariantCaller``

Documentation
-------------


URL
******
*No URL to the documentation was provided*

Docstring
*********
*No documentation was provided: `contribute one <https://github.com/illusional>`_*

Outputs
-------
========  ==========  ===============
name      type        documentation
========  ==========  ===============
splitOut  VCF
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

=======  =======================  ========  ==========  ===============
name     type                     prefix    position    documentation
=======  =======================  ========  ==========  ===============
filters  Optional<Array<String>>
=======  =======================  ========  ==========  ===============


Metadata
********

Author: **Unknown**


*Strelka Variant Caller was last updated on **Unknown***.
*This page was automatically generated on 2019-05-03*.
