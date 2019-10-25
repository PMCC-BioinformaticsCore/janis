:orphan:

Vardict Somatic Variant Caller
============================================================

0 contributors · 1 version

:ID: ``vardictSomaticVariantCaller``
:Python: ``janis_bioinformatics.tools.variantcallers.vardictsomatic_variants import VardictSomaticVariantCaller``
:Versions: v0.1.0
:Authors: 
:Citations: 
:Created: None
:Updated: None
:Required inputs:
   - ``normalBam: BamPair``

   - ``tumorBam: BamPair``

   - ``normalName: String``

   - ``tumorName: String``

   - ``intervals: bed``

   - ``headerLines: File``

   - ``reference: FastaWithDict``
:Outputs: 
   - ``vardictVariants: CompressedVCF``

   - ``out: VCF``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Embedded Tools
***************

======================  ============================
Vardict (Somatic)       ``vardict_somatic/1.6.0``
BCFTools: Annotate      ``bcftoolsAnnotate/v1.5``
Split Multiple Alleles  ``SplitMultiAllele/v0.5772``
Trim IUPAC Bases        ``trimIUPAC/0.0.4``
======================  ============================

------

Additional configuration (inputs)
---------------------------------

============================  =================  ===============
name                          type               documentation
============================  =================  ===============
normalBam                     BamPair
tumorBam                      BamPair
normalName                    String
tumorName                     String
intervals                     bed
headerLines                   File
reference                     FastaWithDict
alleleFreqThreshold           Optional<Float>
vardict_chromNamesAreNumbers  Optional<Boolean>
vardict_vcfFormat             Optional<Boolean>
vardict_chromColumn           Optional<Integer>
vardict_regStartCol           Optional<Integer>
vardict_geneEndCol            Optional<Integer>
============================  =================  ===============

.
