:orphan:

VcfLib: VcfUniqAlleles
=======================================

*1 contributor · 1 version*

:ID: ``vcfuniqalleles``
:Python: ``janis_bioinformatics.tools.vcflib.vcfuniqalleles.versions import VcfUniqAlleles_1_0_1``
:Versions: v1.0.1
:Container: shollizeck/vcflib:1.0.1
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2019-10-18
:Updated: 2019-10-18
:Required inputs:
   - ``vcf: CompressedVCF``
:Outputs: 
   - ``out: stdout<VCF>``

Documentation
-------------

URL: `https://github.com/vcflib/vcflib <https://github.com/vcflib/vcflib>`_

usage: vcffuniq [file]
For each record, remove any duplicate alternate alleles that may have resulted from merging separate VCF files.

------

None

Additional configuration (inputs)
---------------------------------

======  =============  ========  ==========  ===============
name    type           prefix      position  documentation
======  =============  ========  ==========  ===============
vcf     CompressedVCF                     3
======  =============  ========  ==========  ===============

