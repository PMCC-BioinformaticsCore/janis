:orphan:

WGS Somatic (Multi callers)
====================================================

A somatic tumor-normal variant-calling WGS pipeline using GATK, VarDict and Strelka2 · 3 contributors · 1 version

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



4. Generate an inputs file for WGSSomaticMultiCallers:

.. code-block:: bash
   
   janis inputs WGSSomaticMultiCallers > inputs.yaml

**inputs.yaml**

.. code-block:: yaml

       gatk_intervals:
       - gatk_intervals_0.bed
       - gatk_intervals_1.bed
       gridss_blacklist: gridss_blacklist.bed
       known_indels: Homo_sapiens_assembly38.known_indels.vcf.gz
       mills_indels: Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
       normal_inputs:
       - - normal_sample1_R1.fastq.gz
         - normal_sample1_R2.fastq.gz
       - - normal_sample1_R1-TOPUP.fastq.gz
         - normal_sample1_R2-TOPUP.fastq.gz
       normal_name: <value>
       reference: Homo_sapiens_assembly38.fasta
       snps_1000gp: 1000G_phase1.snps.high_confidence.hg38.vcf.gz
       snps_dbsnp: Homo_sapiens_assembly38.dbsnp138.vcf.gz
       strelka_intervals: strelka_intervals.bed.gz
       tumor_inputs:
       - - tumor_sample1_R1.fastq.gz
         - tumor_sample1_R2.fastq.gz
       - - tumor_sample1_R1-TOPUP.fastq.gz
         - tumor_sample1_R2-TOPUP.fastq.gz
       tumor_name: <value>
       vardict_intervals:
       - vardict_intervals_0.bed
       - vardict_intervals_1.bed


5. Run the WGSSomaticMultiCallers pipeline with:

.. code-block:: bash

   janis run [...workflow options] --inputs inputs.yaml WGSSomaticMultiCallers



Outputs
-----------

=================  =================  ================================================
name               type               documentation
=================  =================  ================================================
normal_report      Array<Array<Zip>>  A zip file of the NORMAL FastQC quality reports.
tumor_report       Array<Array<Zip>>  A zip file of the TUMOR FastQC quality reports.
normal_bam         IndexedBam         Aligned and indexed NORMAL bam
tumor_bam          IndexedBam         Aligned and indexed TUMOR bam
gridss_assembly    VCF                Assembly returned by GRIDSS
variants_gatk      VCF                Merged variants from the GATK caller
variants_strelka   VCF                Variants from the Strelka variant caller
variants_vardict   VCF                Merged variants from the VarDict caller
variants_gridss    VCF                Variants from the GRIDSS variant caller
variants_combined  VCF                Combined variants from all 3 callers
=================  =================  ================================================


Information
------------

:ID: ``WGSSomaticMultiCallers``
:Python: ``janis_pipelines.wgs_somatic.wgssomatic import WGSSomaticMultiCallers``
:Versions: 1.2.0
:Authors: Michael Franklin, Richard Lupat, Jiaan Yu
:Citations: 
:Created: 2018-12-24
:Updated: 2020-03-05

Embedded Tools
~~~~~~~~~~~~~~~~~

==============================  ======================================================================================================================================
                                ``somatic_subpipeline/<bound method WorkflowBuilder.version of <janis_core.workflow.workflow.WorkflowBuilder object at 0x1084dccc0>>``
                                ``somatic_subpipeline/<bound method WorkflowBuilder.version of <janis_core.workflow.workflow.WorkflowBuilder object at 0x1084f0320>>``
GATK4 Somatic Variant Caller    ``GATK4_SomaticVariantCaller/4.1.3.0``
GATK4: Gather VCFs              ``Gatk4GatherVcfs/4.1.3.0``
Strelka Somatic Variant Caller  ``strelkaSomaticVariantCaller/v0.1.0``
Gridss                          ``gridss/v2.5.1-dev``
Vardict Somatic Variant Caller  ``vardictSomaticVariantCaller/v0.1.0``
Combine Variants                ``combinevariants/0.0.4``
BCFTools: Sort                  ``bcftoolssort/v1.9``
==============================  ======================================================================================================================================


Additional configuration (inputs)
---------------------------------

========================  =======================  =======================================================================================================================================================================================================================================================================
name                      type                     documentation
========================  =======================  =======================================================================================================================================================================================================================================================================
normal_inputs             Array<FastqGzPair>       An array of NORMAL FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads
tumor_inputs              Array<FastqGzPair>       An array of TUMOR FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads
normal_name               String                   Sample name for the NORMAL sample from which to generate the readGroupHeaderLine for BwaMem
tumor_name                String                   Sample name for the TUMOR sample from which to generate the readGroupHeaderLine for BwaMem
gridss_blacklist          bed
gatk_intervals            Array<bed>               List of intervals over which to split the GATK variant calling
vardict_intervals         Array<bed>               List of intervals over which to split the VarDict variant calling
strelka_intervals         BedTABIX                 An interval for which to restrict the analysis to. Recommended HG38 interval: TBA
reference                 FastaWithIndexes         The reference genome from which to align the reads. This requires a number indexes (can be generated with the 'IndexFasta' pipeline. This pipeline has been tested with the hg38 reference genome.
snps_dbsnp                CompressedIndexedVCF     From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
snps_1000gp               CompressedIndexedVCF     From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
known_indels              CompressedIndexedVCF     From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
mills_indels              CompressedIndexedVCF     From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
cutadapt_adapters         Optional<File>           Specifies a file which contains a list of sequences to determine valid overrepresented sequences from the FastQC report to trim with Cuatadapt. The file must contain sets of named adapters in the form name[tab]sequence. Lines prefixed with a hash will be ignored.
header_lines              Optional<File>           Header lines passed to BCFTools annotate as ``--header-lines``.
allele_freq_threshold     Optional<Float>          The threshold for VarDict's allele frequency, default: 0.05 or 5%
combine_variants_type     Optional<String>         germline | somatic
combine_variants_columns  Optional<Array<String>>  Columns to keep, seperated by space output vcf (unsorted)
========================  =======================  =======================================================================================================================================================================================================================================================================
