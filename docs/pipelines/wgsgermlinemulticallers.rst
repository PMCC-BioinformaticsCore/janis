:orphan:

WGS Germline (Multi callers)
======================================================

A variant-calling WGS pipeline using GATK, VarDict and Strelka2 · 1 contributor · 1 version

:ID: ``WGSGermlineMultiCallers``
:Python: ``janis_pipelines.wgs_germline.wgsgermline import WGSGermlineMultiCallers``
:Versions: 1.2.0
:Authors: Michael Franklin
:Citations: 
:Created: None
:Updated: 2019-10-16
:Required inputs:
   - ``fastqs: Array<FastqGzPair>``

   - ``reference: FastaWithDict``

   - ``cutadapt_adapters: File``

   - ``gatk_intervals: Array<bed>``

   - ``vardict_intervals: Array<bed>``

   - ``strelkaIntervals: BedTABIX``

   - ``header_lines: File``

   - ``snps_dbsnp: CompressedIndexedVCF``

   - ``snps_1000gp: CompressedIndexedVCF``

   - ``known_indels: CompressedIndexedVCF``

   - ``mills_indels: CompressedIndexedVCF``
:Outputs: 
   - ``reports: Array<Array<Zip>>``

   - ``bam: BamPair``

   - ``variants_combined: CompressedVCF``

   - ``variants_gatk: VCF``

   - ``variants_vardict: VCF``

   - ``variants_strelka: VCF``

   - ``variants_gatk_split: Array<VCF>``

   - ``variants_vardict_split: Array<VCF>``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Embedded Tools
***************

===============================  =======================================
FastQC                           ``fastqc/v0.11.5``
Parse FastQC Adaptors            ``ParseFastqcAdaptors/v0.1.0``
Align and sort reads             ``BwaAligner/1.0.0``
Merge and Mark Duplicates        ``mergeAndMarkBams/4.1.3``
GATK4 Germline Variant Caller    ``GATK4_GermlineVariantCaller/4.1.3.0``
GATK4: Gather VCFs               ``Gatk4GatherVcfs/4.1.3.0``
Strelka Germline Variant Caller  ``strelkaGermlineVariantCaller/v0.1.0``
Vardict Germline Variant Caller  ``vardictGermlineVariantCaller/v0.1.0``
Combine Variants                 ``combinevariants/0.0.4``
BCFTools: Sort                   ``bcftoolssort/v1.9``
===============================  =======================================

------

Additional configuration (inputs)
---------------------------------

=============================  =======================  ===============
name                           type                     documentation
=============================  =======================  ===============
fastqs                         Array<FastqGzPair>
reference                      FastaWithDict
cutadapt_adapters              File
gatk_intervals                 Array<bed>
vardict_intervals              Array<bed>
strelkaIntervals               BedTABIX
header_lines                   File
snps_dbsnp                     CompressedIndexedVCF
snps_1000gp                    CompressedIndexedVCF
known_indels                   CompressedIndexedVCF
mills_indels                   CompressedIndexedVCF
sample_name                    Optional<String>
allele_freq_threshold          Optional<Float>
align_and_sort_sortsam_tmpDir  Optional<String>
combine_variants_type          Optional<String>
combine_variants_columns       Optional<Array<String>>
=============================  =======================  ===============

.
