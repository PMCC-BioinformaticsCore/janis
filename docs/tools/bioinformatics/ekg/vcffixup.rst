:orphan:

VcfLib: VcfFixUp
===========================

*1 contributor · 1 version*

:ID: ``vcffixup``
:Python: ``janis_bioinformatics.tools.vcflib.vcffixup.versions import VcfFixUp_1_0_1``
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

usage: vcffixup [file]
Count the allele frequencies across alleles
 present in each record in the VCF file. (Similar to vcftools --freq.)

Uses genotypes from the VCF file to correct AC (alternate allele count), AF (alternate allele frequency), NS (number of called), in the VCF records.

------

None

Additional configuration (inputs)
---------------------------------

======  =============  ========  ==========  ===============
name    type           prefix      position  documentation
======  =============  ========  ==========  ===============
vcf     CompressedVCF                     3
======  =============  ========  ==========  ===============

