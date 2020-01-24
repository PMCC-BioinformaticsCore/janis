:orphan:

Strelka Germline Variant Caller
==============================================================

0 contributors · 1 version

:ID: ``strelkaGermlineVariantCaller``
:Python: ``janis_bioinformatics.tools.variantcallers.illuminagermline_strelka import IlluminaGermlineVariantCaller``
:Versions: v0.1.0
:Authors: 
:Citations: 
:Created: None
:Updated: None
:Required inputs:
   - ``bam: BamPair``

   - ``reference: FastaWithDict``
:Outputs: 
   - ``diploid: CompressedIndexedVCF``

   - ``variants: CompressedIndexedVCF``

   - ``out: VCF``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Embedded Tools
***************

======================  ============================
Manta                   ``manta/1.5.0``
Strelka (Germline)      ``strelka_germline/2.9.10``
BCFTools: View          ``bcftoolsview/v1.5``
Split Multiple Alleles  ``SplitMultiAllele/v0.5772``
======================  ============================

------

Additional configuration (inputs)
---------------------------------

====================  =======================  ===============
name                  type                     documentation
====================  =======================  ===============
bam                   BamPair
reference             FastaWithDict
intervals             Optional<BedTABIX>
is_exome              Optional<Boolean>
bcfview_applyFilters  Optional<Array<String>>
====================  =======================  ===============

.
