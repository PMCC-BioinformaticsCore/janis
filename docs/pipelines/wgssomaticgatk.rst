:orphan:

WGS Somatic (GATK only)
========================================

*A somatic tumor-normal variant-calling WGS pipeline using only GATK Mutect2 · 3 contributors · 1 version*

This is a genomics pipeline to align sequencing data (Fastq pairs) into BAMs:

- Takes raw sequence data in the FASTQ format;
- align to the reference genome using BWA MEM;
- Marks duplicates using Picard;
- Call the appropriate somatic variant callers (GATK / Strelka / VarDict);
- Outputs the final variants in the VCF format.

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

=================  ====================  =========================================================================================================================================================================================  ======================================================================================================================================================================================================================================================================================================
Name               Type                  Example                                                                                                                                                                                    Description
=================  ====================  =========================================================================================================================================================================================  ======================================================================================================================================================================================================================================================================================================
cutadapt_adapters  Optional<File>        https://github.com/csf-ngs/fastqc/blob/master/Contaminants/contaminant_list.txt                                                                                                            Specifies a containment list for cutadapt, which contains a list of sequences to determine valid overrepresented sequences from the FastQC report to trim with Cuatadapt. The file must contain sets of named adapters in the form: ``name[tab]sequence``. Lines prefixed with a hash will be ignored.
gatk_intervals     Array<bed>            BRCA1.bed                                                                                                                                                                                  List of intervals over which to split the GATK variant calling
reference          FastaWithIndexes      HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                                                                            The reference genome from which to align the reads. This requires a number indexes (can be generated with the 'IndexFasta' pipeline This pipeline has been tested using the HG38 reference set.

                                         File: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta                                                                                                           This pipeline expects the assembly references to be as they appear in the GCP example:

                                                                                                                                                                                                                                    - (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
snps_dbsnp         CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                                                                            From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                         (WARNING: The file available from the genomics-public-data resource on Google Cloud Storage is NOT compressed and indexed. This will need to be completed prior to starting the pipeline.

                                         File: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
snps_1000gp        CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                                                                            From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                         File: gs://genomics-public-data/references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
known_indels       CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                                                                            From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                         File: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
mills_indels       CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                                                                            From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                         File: gs://genomics-public-data/references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
=================  ====================  =========================================================================================================================================================================================  ======================================================================================================================================================================================================================================================================================================

4. Generate user and static input files for WGSSomaticGATK:

.. code-block:: bash

   # user inputs
   janis inputs --user WGSSomaticGATK > inputs.yaml

   # static inputs
   janis inputs --static WGSSomaticGATK > static.yaml

**inputs.yaml**

.. code-block:: yaml

       normal_inputs:
       - - normal_R1.fastq.gz
         - normal_R2.fastq.gz
       - - normal_R1-TOPUP.fastq.gz
         - normal_R2-TOPUP.fastq.gz
       normal_name: <value>
       tumor_inputs:
       - - tumor_R1.fastq.gz
         - tumor_R2.fastq.gz
       - - tumor_R1-TOPUP.fastq.gz
         - tumor_R2-TOPUP.fastq.gz
       tumor_name: <value>


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


5. Run WGSSomaticGATK with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       --inputs static.yaml \
       WGSSomaticGATK



Outputs
-----------

=============  =================  ===============
name           type               documentation
=============  =================  ===============
normal_bam     IndexedBam
tumor_bam      IndexedBam
normal_report  Array<Array<Zip>>
tumor_report   Array<Array<Zip>>
variants_gatk  CompressedVCF
=============  =================  ===============


Information
------------

:ID: ``WGSSomaticGATK``
:Versions: 1.2.0
:Authors: Michael Franklin, Richard Lupat, Jiaan Yu
:Citations: 
:Created: None
:Updated: 2020-03-16

Embedded Tools
~~~~~~~~~~~~~~~~~

============================  ======================================================================================================================================
                              ``somatic_subpipeline/<bound method WorkflowBuilder.version of <janis_core.workflow.workflow.WorkflowBuilder object at 0x11258d860>>``
                              ``somatic_subpipeline/<bound method WorkflowBuilder.version of <janis_core.workflow.workflow.WorkflowBuilder object at 0x11259db38>>``
GATK4 Somatic Variant Caller  ``GATK4_SomaticVariantCaller/4.1.3.0``
GATK4: Gather VCFs            ``Gatk4GatherVcfs/4.1.3.0``
BCFTools: Sort                ``bcftoolssort/v1.9``
============================  ======================================================================================================================================


Additional configuration (inputs)
---------------------------------

=================  ====================  ======================================================================================================================================================================================================================================================================================================
name               type                  documentation
=================  ====================  ======================================================================================================================================================================================================================================================================================================
normal_inputs      Array<FastqGzPair>    An array of NORMAL FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads
tumor_inputs       Array<FastqGzPair>    An array of TUMOR FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads
normal_name        String                Sample name for the NORMAL sample from which to generate the readGroupHeaderLine for BwaMem
tumor_name         String                Sample name for the TUMOR sample from which to generate the readGroupHeaderLine for BwaMem
gatk_intervals     Array<bed>            List of intervals over which to split the GATK variant calling
reference          FastaWithIndexes      The reference genome from which to align the reads. This requires a number indexes (can be generated with the 'IndexFasta' pipeline This pipeline has been tested using the HG38 reference set.

                                         This pipeline expects the assembly references to be as they appear in the GCP example:

                                         - (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
snps_dbsnp         CompressedIndexedVCF  From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
snps_1000gp        CompressedIndexedVCF  From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
known_indels       CompressedIndexedVCF  From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
mills_indels       CompressedIndexedVCF  From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
cutadapt_adapters  Optional<File>        Specifies a containment list for cutadapt, which contains a list of sequences to determine valid overrepresented sequences from the FastQC report to trim with Cuatadapt. The file must contain sets of named adapters in the form: ``name[tab]sequence``. Lines prefixed with a hash will be ignored.
=================  ====================  ======================================================================================================================================================================================================================================================================================================
