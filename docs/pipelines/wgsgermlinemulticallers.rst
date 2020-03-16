:orphan:

WGS Germline (Multi callers)
======================================================

*A variant-calling WGS pipeline using GATK, VarDict and Strelka2 · 3 contributors · 1 version*

This is a genomics pipeline to align sequencing data (Fastq pairs) into BAMs and call variants using:

This workflow is a reference pipeline using the Janis Python framework (pipelines assistant).

- Takes raw sequence data in the FASTQ format;
- align to the reference genome using BWA MEM;
- Marks duplicates using Picard;
- Call the appropriate variant callers (GATK / Strelka / VarDict);
- Outputs the final variants in the VCF format.


Quickstart
-----------

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.

====================  ====================  ===========================================================================================================================================  =================================================================================================================================================================================================================================================================================================================================================================================
Name                  Type                  Example                                                                                                                                      Description
====================  ====================  ===========================================================================================================================================  =================================================================================================================================================================================================================================================================================================================================================================================
reference             FastaWithIndexes      HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                              The reference genome from which to align the reads. This requires a number indexes (can be generated with the 'IndexFasta' pipeline This pipeline has been tested using the HG38 reference set.

                                            File: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta                                                             This pipeline expects the assembly references to be as they appear in the GCP example:

                                                                                                                                                                                         - (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
cutadapt_adapters     Optional<File>        https://github.com/csf-ngs/fastqc/blob/master/Contaminants/contaminant_list.txt                                                              Specifies a containment list for cutadapt, which contains a list of sequences to determine valid overrepresented sequences from the FastQC report to trim with Cuatadapt. The file must contain sets of named adapters in the form: ``name[tab]sequence``. Lines prefixed with a hash will be ignored.
gatk_intervals        Array<bed>            BRCA1.bed                                                                                                                                    List of intervals over which to split the GATK variant calling
vardict_intervals     Array<bed>            BRCA1.bed                                                                                                                                    List of intervals over which to split the VarDict variant calling
strelka_intervals     BedTABIX              BRCA1.bed.gz                                                                                                                                 An interval for which to restrict the analysis to.
vardict_header_lines  File                  https://gist.githubusercontent.com/illusional/5b75a0506f7327aca7d355f8ad5008f8/raw/e181c0569771e6a557d01a8a1f70c71e3598a269/headerLines.txt  As with chromosomal sequences it is highly recommended (but not required) that the header include tags describing the contigs referred to in the VCF file. This furthermore allows these contigs to come from different files. The format is identical to that of a reference sequence, but with an additional URL tag to indicate where that sequence can be found. For example:

                                                                                                                                                                                         .. code-block:

                                                                                                                                                                                            ##contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>

                                                                                                                                                                                         Source: (1.2.5 Alternative allele field format) https://samtools.github.io/hts-specs/VCFv4.1.pdf (edited)
snps_dbsnp            CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                              From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                            File: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
snps_1000gp           CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                              From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                            File: gs://genomics-public-data/references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
known_indels          CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                              From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                            File: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
mills_indels          CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                              From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                            File: gs://genomics-public-data/references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
====================  ====================  ===========================================================================================================================================  =================================================================================================================================================================================================================================================================================================================================================================================

4. Generate user and static input files for WGSGermlineMultiCallers:

.. code-block:: bash
   
   # user inputs
   janis inputs --user WGSGermlineMultiCallers > inputs.yaml
    
   # static inputs
   janis inputs --static WGSGermlineMultiCallers > static.yaml


**inputs.yaml**

.. code-block:: yaml

       fastqs:
       - - sample1_R1.fastq.gz
         - sample1_R2.fastq.gz
       - - sample1_R1-TOPUP.fastq.gz
         - sample1_R2-TOPUP.fastq.gz
       sample_name: <value>


**static.yaml**

.. code-block:: yaml

       gatk_intervals:
       - gatk_intervals_0.bed
       - gatk_intervals_1.bed
       known_indels: Homo_sapiens_assembly38.known_indels.vcf.gz
       mills_indels: Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
       reference: Homo_sapiens_assembly38.fasta
       snps_1000gp: 1000G_phase1.snps.high_confidence.hg38.vcf.gz
       snps_dbsnp: Homo_sapiens_assembly38.dbsnp138.vcf.gz
       strelka_intervals: strelka_intervals.bed.gz
       vardict_header_lines: vardict_header_lines
       vardict_intervals:
       - vardict_intervals_0.bed
       - vardict_intervals_1.bed


5. Run the WGSGermlineMultiCallers pipeline with:

.. code-block:: bash

   janis run [...workflow options] \
       --inputs inputs.yaml \
       --inputs static.yaml \
       WGSGermlineMultiCallers



Outputs
-----------

======================  =================  =======================================================
name                    type               documentation
======================  =================  =======================================================
reports                 Array<Array<Zip>>  A zip file of the FastQC quality report.
bam                     IndexedBam         Aligned and indexed bam.
variants_combined       CompressedVCF      Combined variants from all 3 callers
variants_gatk           VCF                Merged variants from the GATK caller
variants_vardict        VCF                Merged variants from the VarDict caller
variants_strelka        VCF                Variants from the Strelka variant caller
variants_gatk_split     Array<VCF>         Unmerged variants from the GATK caller (by interval)
variants_vardict_split  Array<VCF>         Unmerged variants from the VarDict caller (by interval)
======================  =================  =======================================================


Information
------------

:ID: ``WGSGermlineMultiCallers``
:Python: ``janis_pipelines.wgs_germline.wgsgermline import WGSGermlineMultiCallers``
:Versions: 1.2.0
:Authors: Michael Franklin, Richard Lupat, Jiaan Yu
:Citations: 
:Created: 2018-12-24
:Updated: 2020-03-16

Embedded Tools
~~~~~~~~~~~~~~~~~

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


Additional configuration (inputs)
---------------------------------

=============================  =======================  =================================================================================================================================================================================================================================================================================================================================================================================
name                           type                     documentation
=============================  =======================  =================================================================================================================================================================================================================================================================================================================================================================================
sample_name                    String                   Sample name from which to generate the readGroupHeaderLine for BwaMem
fastqs                         Array<FastqGzPair>       An array of FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads
reference                      FastaWithIndexes         The reference genome from which to align the reads. This requires a number indexes (can be generated with the 'IndexFasta' pipeline This pipeline has been tested using the HG38 reference set.

                                                        This pipeline expects the assembly references to be as they appear in the GCP example:

                                                        - (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
gatk_intervals                 Array<bed>               List of intervals over which to split the GATK variant calling
vardict_intervals              Array<bed>               List of intervals over which to split the VarDict variant calling
strelka_intervals              BedTABIX                 An interval for which to restrict the analysis to.
vardict_header_lines           File                     As with chromosomal sequences it is highly recommended (but not required) that the header include tags describing the contigs referred to in the VCF file. This furthermore allows these contigs to come from different files. The format is identical to that of a reference sequence, but with an additional URL tag to indicate where that sequence can be found. For example:

                                                        .. code-block:

                                                           ##contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>

                                                        Source: (1.2.5 Alternative allele field format) https://samtools.github.io/hts-specs/VCFv4.1.pdf (edited)
snps_dbsnp                     CompressedIndexedVCF     From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
snps_1000gp                    CompressedIndexedVCF     From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
known_indels                   CompressedIndexedVCF     From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
mills_indels                   CompressedIndexedVCF     From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
cutadapt_adapters              Optional<File>           Specifies a containment list for cutadapt, which contains a list of sequences to determine valid overrepresented sequences from the FastQC report to trim with Cuatadapt. The file must contain sets of named adapters in the form: ``name[tab]sequence``. Lines prefixed with a hash will be ignored.
allele_freq_threshold          Optional<Float>          The threshold for VarDict's allele frequency, default: 0.05 or 5%
align_and_sort_sortsam_tmpDir  Optional<String>         Undocumented option
combine_variants_type          Optional<String>         germline | somatic
combine_variants_columns       Optional<Array<String>>  Columns to keep, seperated by space output vcf (unsorted)
=============================  =======================  =================================================================================================================================================================================================================================================================================================================================================================================
