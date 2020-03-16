:orphan:

Vardict Germline Variant Caller
==============================================================

*0 contributors · 1 version*

:ID: ``vardictGermlineVariantCaller``
:Python: ``janis_bioinformatics.tools.variantcallers.vardictgermline_variants import VardictGermlineVariantCaller``
:Versions: v0.1.0
:Authors: 
:Citations: 
:Created: None
:Updated: None
:Required inputs:
   - ``bam: IndexedBam``

   - ``intervals: bed``

   - ``sample_name: String``

   - ``header_lines: File``

   - ``reference: FastaWithIndexes``
:Outputs: 
   - ``vardict_variants: CompressedVCF``

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
Trim IUPAC Bases        ``trimIUPAC/0.0.5``
======================  ============================

------

Additional configuration (inputs)
---------------------------------

============================  =================  ============================================================================
name                          type               documentation
============================  =================  ============================================================================
bam                           IndexedBam
intervals                     bed
sample_name                   String
header_lines                  File
reference                     FastaWithIndexes
allele_freq_threshold         Optional<Float>
vardict_chromNamesAreNumbers  Optional<Boolean>  Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
vardict_vcfFormat             Optional<Boolean>  VCF format output
vardict_chromColumn           Optional<Integer>  The column for chromosome
vardict_regStartCol           Optional<Integer>  The column for region start, e.g. gene start
vardict_geneEndCol            Optional<Integer>  The column for region end, e.g. gene end
============================  =================  ============================================================================

.
