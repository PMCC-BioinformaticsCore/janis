:orphan:

WGS Germline (GATK only)
==========================================

A variant-calling WGS pipeline using only the GATK Haplotype variant caller · 1 contributor · 1 version

:ID: ``WGSGermlineGATK``
:Python: ``janis_pipelines.wgs_germline_gatk.wgsgermlinegatk import WGSGermlineGATK``
:Versions: 1.0.0
:Contributors: Michael Franklin
:Citations: 
:Created: None
:Updated: 2019-10-16
:Required inputs:
   - ``fastqs: Array<FastqGzPair>``

   - ``reference: FastaWithDict``

   - ``gatkIntervals: Array<bed>``

   - ``snps_dbsnp: CompressedIndexedVCF``

   - ``snps_1000gp: CompressedIndexedVCF``

   - ``known_indels: CompressedIndexedVCF``

   - ``mills_indels: CompressedIndexedVCF``
:Outputs: 
   - ``bam: BamPair``

   - ``reports: Array<Array<Zip>>``

   - ``variants: CompressedVCF``

   - ``variants_split: Array<VCF>``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Embedded Tools
***************

=============================  ======================================
FastQC                         ``fastqc/v0.11.5``
Align and sort reads           ``BwaAligner/1.0.0``
Merge and Mark Duplicates      ``mergeAndMarkBams/4.1.3``
GATK4 Germline Variant Caller  ``GATK4_GermlineVariantCaller/v0.1.0``
GATK4: Gather VCFs             ``Gatk4GatherVcfs/4.0.12.0``
BCFTools: Sort                 ``bcftoolssort/v1.9``
=============================  ======================================

------

Inputs
------

=============================  ====================  ===============
name                           type                  documentation
=============================  ====================  ===============
fastqs                         Array<FastqGzPair>
reference                      FastaWithDict
gatkIntervals                  Array<bed>
snps_dbsnp                     CompressedIndexedVCF
snps_1000gp                    CompressedIndexedVCF
known_indels                   CompressedIndexedVCF
mills_indels                   CompressedIndexedVCF
sampleName                     Optional<String>
alignSortedBam_sortsam_tmpDir  Optional<String>
=============================  ====================  ===============

.
