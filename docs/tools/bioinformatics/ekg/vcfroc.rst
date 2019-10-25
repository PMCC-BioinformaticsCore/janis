:orphan:

VcfLib: Vcf ROC generator
==================================

1 contributor · 1 version

:ID: ``vcfroc``
:Python: ``janis_bioinformatics.tools.vcflib.vcfroc.versions import VcfRoc_1_0_1``
:Versions: v1.0.1
:Container: shollizeck/vcflib:1.0.1
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2019-10-18
:Updated: 2019-10-18
:Required inputs:
   - ``vcf: CompressedVCF``

   - ``truth: CompressedVCF``

   - ``reference: FastaWithDict``
:Outputs: 
   - ``out: stdout<VCF>``

Documentation
-------------

URL: `https://github.com/vcflib/vcflib <https://github.com/vcflib/vcflib>`_

usage: vcfroc [options] [<vcf file>]

options:
	-t, --truth-vcf FILE	use this VCF as ground truth for ROC generation
	-w, --window-size N       compare records up to this many bp away (default 30)
	-r, --reference FILE	FASTA reference file

------

Additional configuration (inputs)
---------------------------------

==========  =================  ====================================================
name        type               documentation
==========  =================  ====================================================
vcf         CompressedVCF
truth       CompressedVCF      use this VCF as ground truth for ROC generation
reference   FastaWithDict      FASTA reference file
windowSize  Optional<Integer>  compare records up to this many bp away (default 30)
==========  =================  ====================================================

