:orphan:

VcfLib: Vcf Random Sampling
=============================================

1 contributor · 1 version

:ID: ``vcfrandomsample``
:Python: ``janis_bioinformatics.tools.vcflib.vcfrandomsample.versions import VcfRandomSample_1_0_1``
:Versions: v1.0.1
:Container: shollizeck/vcflib:1.0.1
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2019-10-18
:Updated: 2019-10-18
:Required inputs:
   - ``vcf: CompressedVCF``

   - ``rate: Float``

   - ``seed: Integer``
:Outputs: 
   - ``out: stdout<VCF>``

Documentation
-------------

URL: `https://github.com/vcflib/vcflib <https://github.com/vcflib/vcflib>`_

usage: vcfrandomsample [options] [<vcf file>]

options:
	-r, --rate RATE 	base sampling probability per locus
	-s, --scale-by KEY\scale sampling likelihood by this Float info field
	-p, --random-seed N	use this random seed

------

Additional configuration (inputs)
---------------------------------

=======  ================  ==================================================
name     type              documentation
=======  ================  ==================================================
vcf      CompressedVCF
rate     Float             base sampling probability per locus
seed     Integer           use this random seed
scaleBy  Optional<String>  scale sampling likelihood by this Float info field
=======  ================  ==================================================

