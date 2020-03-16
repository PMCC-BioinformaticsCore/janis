:orphan:

WGS Germline (GATK only)
==========================================

*A variant-calling WGS pipeline using only the GATK Haplotype variant caller · 3 contributors · 1 version*

This is a genomics pipeline to align sequencing data (Fastq pairs) into BAMs and call variants using GATK. The final variants are outputted in the VCF format.

**Resources**

This pipeline has been tested using the HG38 reference set, available on Google Cloud Storage through:

- https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/

This pipeline expects the assembly references to be as they appear in that storage     (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
The known sites (snps_dbsnp, snps_1000gp, known_indels, mills_indels) should be gzipped and tabix indexed.


Quickstart
-----------

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.

=================  ====================  ===================================================================================================  ======================================================================================================================================================================================================================================================================================================
Name               Type                  Example                                                                                              Description
=================  ====================  ===================================================================================================  ======================================================================================================================================================================================================================================================================================================
reference          FastaWithIndexes      HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/      The reference genome from which to align the reads. This requires a number indexes (can be generated with the 'IndexFasta' pipeline This pipeline has been tested using the HG38 reference set.

                                         File: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta                     This pipeline expects the assembly references to be as they appear in the GCP example:

                                                                                                                                              - (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
cutadapt_adapters  Optional<File>        https://github.com/csf-ngs/fastqc/blob/master/Contaminants/contaminant_list.txt                      Specifies a containment list for cutadapt, which contains a list of sequences to determine valid overrepresented sequences from the FastQC report to trim with Cuatadapt. The file must contain sets of named adapters in the form: ``name[tab]sequence``. Lines prefixed with a hash will be ignored.
gatk_intervals     Array<bed>            BRCA1.bed                                                                                            List of intervals over which to split the GATK variant calling
snps_dbsnp         CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/      From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                         File: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
snps_1000gp        CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/      From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                         File: gs://genomics-public-data/references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
known_indels       CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/      From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                         File: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
mills_indels       CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/      From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                         File: gs://genomics-public-data/references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
=================  ====================  ===================================================================================================  ======================================================================================================================================================================================================================================================================================================

4. Generate user and static input files for WGSGermlineGATK:

.. code-block:: bash

   # user inputs
   janis inputs --user WGSGermlineGATK > inputs.yaml

   # static inputs
   janis inputs --static WGSGermlineGATK > static.yaml

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


5. Run WGSGermlineGATK with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       --inputs static.yaml \
       WGSGermlineGATK



Outputs
-----------

==============  =================  ====================================================
name            type               documentation
==============  =================  ====================================================
bam             IndexedBam         Aligned and indexed bam.
reports         Array<Array<Zip>>  A zip file of the FastQC quality report.
variants        CompressedVCF      Merged variants from the GATK caller
variants_split  Array<VCF>         Unmerged variants from the GATK caller (by interval)
==============  =================  ====================================================


Information
------------

:ID: ``WGSGermlineGATK``
:Versions: 1.2.0
:Authors: Michael Franklin, Richard Lupat, Jiaan Yu
:Citations: 
:Created: 2018-12-24
:Updated: 2020-03-16

Embedded Tools
~~~~~~~~~~~~~~~~~

=============================  =======================================
FastQC                         ``fastqc/v0.11.5``
Parse FastQC Adaptors          ``ParseFastqcAdaptors/v0.1.0``
Align and sort reads           ``BwaAligner/1.0.0``
Merge and Mark Duplicates      ``mergeAndMarkBams/4.1.3``
GATK4 Germline Variant Caller  ``GATK4_GermlineVariantCaller/4.1.3.0``
GATK4: Gather VCFs             ``Gatk4GatherVcfs/4.0.12.0``
BCFTools: Sort                 ``bcftoolssort/v1.9``
=============================  =======================================


Additional configuration (inputs)
---------------------------------

=============================  ====================  ======================================================================================================================================================================================================================================================================================================
name                           type                  documentation
=============================  ====================  ======================================================================================================================================================================================================================================================================================================
sample_name                    String                Sample name from which to generate the readGroupHeaderLine for BwaMem
fastqs                         Array<FastqGzPair>    An array of FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads
reference                      FastaWithIndexes      The reference genome from which to align the reads. This requires a number indexes (can be generated with the 'IndexFasta' pipeline This pipeline has been tested using the HG38 reference set.

                                                     This pipeline expects the assembly references to be as they appear in the GCP example:

                                                     - (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
gatk_intervals                 Array<bed>            List of intervals over which to split the GATK variant calling
snps_dbsnp                     CompressedIndexedVCF  From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
snps_1000gp                    CompressedIndexedVCF  From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
known_indels                   CompressedIndexedVCF  From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
mills_indels                   CompressedIndexedVCF  From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
cutadapt_adapters              Optional<File>        Specifies a containment list for cutadapt, which contains a list of sequences to determine valid overrepresented sequences from the FastQC report to trim with Cuatadapt. The file must contain sets of named adapters in the form: ``name[tab]sequence``. Lines prefixed with a hash will be ignored.
align_and_sort_sortsam_tmpDir  Optional<String>      Undocumented option
=============================  ====================  ======================================================================================================================================================================================================================================================================================================
