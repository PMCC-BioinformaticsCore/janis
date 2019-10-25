:orphan:

VcfLib: VcfAllelicPrimitives
===================================================

1 contributor · 1 version

:ID: ``vcfallelicprimitives``
:Python: ``janis_bioinformatics.tools.vcflib.vcfallelicprimitives.versions import VcfAllelicPrimitives_1_0_1``
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

usage: vcfallelicprimitives [options] [file]

options:
	-m, --use-mnps	Retain MNPs as separate events (default: false)
	-t, --tag-parsed FLAG	Tag records which are split apart of a complex allele with this flag

------

Additional configuration (inputs)
---------------------------------

===========  =================  ====================================================================
name         type               documentation
===========  =================  ====================================================================
vcf          CompressedVCF
useMnpsFlag  Optional<Boolean>  Retain MNPs as separate events (default: false)
tagParsed    Optional<String>   Tag records which are split apart of a complex allele with this flag
===========  =================  ====================================================================

