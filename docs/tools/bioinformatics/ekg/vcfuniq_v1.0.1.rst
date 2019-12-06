:orphan:

VcfLib: VcfUniq
=========================

1 contributor · 1 version

:ID: ``vcfuniq``
:Python: ``janis_bioinformatics.tools.vcflib.vcfuniq.versions import VcfUniq_1_0_1``
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
Like GNU uniq, but for VCF records. Remove records which have the same positon, ref, and alt as the previous record.

------

None

Additional configuration (inputs)
---------------------------------

======  =============  ========  ==========  ===============
name    type           prefix      position  documentation
======  =============  ========  ==========  ===============
vcf     CompressedVCF                     3
======  =============  ========  ==========  ===============

