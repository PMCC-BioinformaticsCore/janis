:orphan:

Gridss Germline Variant Caller
============================================================

0 contributors · 1 version

:ID: ``gridssGermlineVariantCaller``
:Python: ``janis_bioinformatics.tools.variantcallers.gridssgermline import GridssGermlineVariantCaller``
:Versions: v0.1.0
:Authors: 
:Citations: 
:Created: None
:Updated: None
:Required inputs:
   - ``bam: BamPair``

   - ``reference: FastaWithDict``

   - ``blacklist: bed``
:Outputs: 
   - ``out: VCF``

   - ``assembly: BAM``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Embedded Tools
***************

==============  ======================
SamTools: View  ``SamToolsView/1.7.0``
Gridss          ``gridss/v2.5.1-dev``
==============  ======================

------

Additional configuration (inputs)
---------------------------------

=========================================  ================  ===============
name                                       type              documentation
=========================================  ================  ===============
bam                                        BamPair
reference                                  FastaWithDict
blacklist                                  bed
samtools_doNotOutputAlignmentsWithBitsSet  Optional<String>
=========================================  ================  ===============

.
