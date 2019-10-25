:orphan:

Strelka Somatic Variant Caller
============================================================

0 contributors · 1 version

:ID: ``strelkaSomaticVariantCaller``
:Python: ``janis_bioinformatics.tools.variantcallers.illuminasomatic_strelka import IlluminaSomaticVariantCaller``
:Versions: v0.1.0
:Authors: 
:Citations: 
:Created: None
:Updated: None
:Required inputs:
   - ``normalBam: BamPair``

   - ``tumorBam: BamPair``

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
Strelka (Somatic)       ``strelka_somatic/2.9.10``
BCFTools: View          ``bcftoolsview/v1.5``
Split Multiple Alleles  ``SplitMultiAllele/v0.5772``
======================  ============================

------

Additional configuration (inputs)
---------------------------------

=====================  =======================  ===============
name                   type                     documentation
=====================  =======================  ===============
normalBam              BamPair
tumorBam               BamPair
reference              FastaWithDict
intervals              Optional<BedTABIX>
isExome                Optional<Boolean>
bcf_view_applyFilters  Optional<Array<String>>
=====================  =======================  ===============

.
