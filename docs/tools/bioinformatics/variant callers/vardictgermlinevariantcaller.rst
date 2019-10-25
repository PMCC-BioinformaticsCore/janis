:orphan:

Vardict Germline Variant Caller
==============================================================

0 contributors · 1 version

:ID: ``vardictGermlineVariantCaller``
:Python: ``janis_bioinformatics.tools.variantcallers.vardictgermline_variants import VardictGermlineVariantCaller``
:Versions: v0.1.0
:Authors: 
:Citations: 
:Created: None
:Updated: None
:Required inputs:
   - ``bam: BamPair``

   - ``intervals: bed``

   - ``sampleName: String``

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
VarDict (Germline)      ``vardict_germline/1.6.0``
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
bam                           BamPair
intervals                     bed
sampleName                    String
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
